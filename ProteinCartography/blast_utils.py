#!/usr/bin/env python
"""BLAST runner for ProteinCartography.

Stock ``blastp -remote`` often hangs from cloud IPs with no timeout. This module
supports a selectable backend:

``PC_BLAST_BACKEND``
  - ``www`` (default) — NCBI Blast.cgi URL API (Put → RID poll → Get XML),
    then write the same tabular TSV ProteinCartography expects
    (``BLAST_OUTFMT`` / ``sacc`` column). This is the systematic fix.
  - ``remote`` — stock ``blastp -remote`` CLI, with a hard timeout + process-group
    kill (last resort / local debugging).
  - ``auto`` — try ``www``, then ``remote``.

Env knobs
---------
``PC_BLAST_TIMEOUT`` — overall wall clock (default 1200s for www).
``PC_BLAST_EMAIL`` / ``PC_BLAST_TOOL`` — NCBI usage-guideline identity.
``PC_BLAST_POLL_SECONDS`` — RID poll interval (NCBI asks ≥60s; default 60).
"""

from __future__ import annotations
import os
import re
import subprocess
import threading
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def blast_call_failed(result) -> bool:
    """
    Determine whether a call to the `blastp` CLI failed.

    The exit code alone is not sufficient: when the remote server refuses to queue a request,
    `blastp -remote` writes an error to stderr, writes an empty results file, and nevertheless
    exits with a status of zero. For example:

        Error: [blastp] bad_request: Could not queue request: DB operation failed.

    Checking only the exit code therefore treats these failures as successes, which prevents the
    word-size backoff in `run_blast.py` from ever being attempted and defers the failure to
    `extract_blast_hits.py`, where it surfaces as a misleading "no hits were returned".

    Note: an empty results file is deliberately *not* treated as a failure here, because a query
    that legitimately has no hits also produces one.

    Args:
        result: the result of the call, with `returncode` and `stderr` attributes.

    Returns:
        True if the call failed.
    """
    if result.returncode != 0:
        return True

    stderr = result.stderr
    if isinstance(stderr, bytes):
        stderr = stderr.decode(errors="replace")

    return "Error:" in (stderr or "")


# Must match ``constants.BLAST_OUTPUT_FIELDS`` order after the leading ``6``
# in ``BLAST_OUTFMT``.
BLAST_TSV_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "sacc",
    "saccver",
    "sgi",
    "staxids",
    "scomnames",
]


@dataclass
class BlastResult:
    returncode: int
    stdout: bytes
    stderr: bytes


def run_blast(
    query: str,
    out: str,
    max_target_seqs: int,
    outfmt: str,
    word_size: int,
    evalue: float,
):
    """
    Call BLAST against NCBI and write tabular results to ``out``.

    Argument names match (a subset of) the blastp CLI / PC ``run_blast.py``.
    Never returns ``None``.
    """
    backend = os.environ.get("PC_BLAST_BACKEND", "www").strip().lower() or "www"
    if backend not in {"www", "remote", "auto"}:
        print(f"[blast] unknown PC_BLAST_BACKEND={backend!r}; using www", flush=True)
        backend = "www"

    result: BlastResult | None = None
    try:
        if backend == "www":
            result = _run_blast_www(query, out, max_target_seqs, word_size, evalue)
        elif backend == "remote":
            result = _run_blast_remote(query, out, max_target_seqs, outfmt, word_size, evalue)
        else:
            result = _run_blast_www(query, out, max_target_seqs, word_size, evalue)
            ok = result is not None and result.returncode == 0 and Path(out).exists()
            if not ok:
                print(
                    "[blast] www backend failed; falling back to blastp -remote",
                    flush=True,
                )
                result = _run_blast_remote(query, out, max_target_seqs, outfmt, word_size, evalue)
    except Exception as exc:
        msg = f"[blast] unexpected error in run_blast: {type(exc).__name__}: {exc}\n"
        print(msg, flush=True)
        return BlastResult(returncode=1, stdout=b"", stderr=msg.encode())

    if result is None:
        msg = "[blast] backend returned None (treating as failure)\n"
        print(msg, flush=True)
        return BlastResult(returncode=1, stdout=b"", stderr=msg.encode())
    return result


# --------------------------------------------------------------------------- #
# NCBI WWW / QBlast (systematic fix)
# --------------------------------------------------------------------------- #


def _run_blast_www(
    query: str,
    out: str,
    max_target_seqs: int,
    word_size: int,
    evalue: float,
) -> BlastResult:
    timeout = float(os.environ.get("PC_BLAST_TIMEOUT", "1200"))
    poll = float(os.environ.get("PC_BLAST_POLL_SECONDS", "60"))
    # NCBI: do not poll a RID more often than once a minute.
    poll = max(60.0, poll)
    email = os.environ.get("PC_BLAST_EMAIL", "proteincartography@arcadia.science").strip()
    tool = os.environ.get("PC_BLAST_TOOL", "ProteinCartography").strip()
    # Default to ``refseq_protein``: NCBI often queues ``nr`` for hours from cloud
    # IPs, and map_refseq_ids already expects RefSeq accessions. Override with
    # ``PC_BLAST_DB=nr`` for the historical broader search.
    default_db = "refseq_protein"
    database = os.environ.get("PC_BLAST_DB", default_db).strip() or default_db

    try:
        fasta = Path(query).read_text()
    except OSError as exc:
        msg = f"[blast] cannot read query FASTA {query}: {exc}\n"
        print(msg, flush=True)
        return BlastResult(returncode=1, stdout=b"", stderr=msg.encode())

    print(
        f"[blast] www/QBlast submit db={database} hitlist={max_target_seqs} "
        f"word_size={word_size} evalue={evalue} timeout={timeout:.0f}s "
        f"query={query}",
        flush=True,
    )

    stop = threading.Event()
    started = time.time()

    def _heartbeat():
        while not stop.wait(30.0):
            elapsed = time.time() - started
            print(
                f"[blast] www still waiting on NCBI ({elapsed:.0f}s / {timeout:.0f}s)",
                flush=True,
            )

    hb = threading.Thread(target=_heartbeat, name="blast-www-hb", daemon=True)
    hb.start()
    try:
        rid = None
        rtoe = 30
        put_errors: list[str] = []
        for put_try in range(1, 4):
            try:
                rid, rtoe = _www_put(
                    fasta,
                    database=database,
                    hitlist_size=max_target_seqs,
                    word_size=word_size,
                    evalue=evalue,
                    email=email,
                    tool=tool,
                )
                break
            except Exception as exc:
                put_errors.append(str(exc))
                print(f"[blast] www PUT try {put_try}/3 failed: {exc}", flush=True)
                time.sleep(10 * put_try)
        if rid is None:
            msg = f"[blast] www PUT failed after retries: {put_errors[-1:]}\n"
            print(msg, flush=True)
            return BlastResult(returncode=1, stdout=b"", stderr=msg.encode())

        print(f"[blast] www RID={rid} RTOE≈{rtoe}s; polling every {poll:.0f}s", flush=True)
        # Honour NCBI's own estimate (capped) so we don't soft-fail while the
        # job is still legitimately queued.
        effective_timeout = max(timeout, min(float(rtoe or 0) + 300.0, 7200.0))
        if effective_timeout > timeout:
            print(
                f"[blast] www extending timeout {timeout:.0f}s → {effective_timeout:.0f}s "
                f"from RTOE",
                flush=True,
            )
        initial_wait = min(max(float(rtoe or 0), 0.0), 120.0)
        if initial_wait > 0:
            time.sleep(initial_wait)

        while True:
            elapsed = time.time() - started
            if elapsed > effective_timeout:
                msg = (
                    f"[blast] www TIMEOUT after {effective_timeout:.0f}s (RID={rid} RTOE≈{rtoe}s)\n"
                )
                print(msg, flush=True)
                return BlastResult(returncode=124, stdout=b"", stderr=msg.encode())
            try:
                status, body = _www_get(rid, email=email, tool=tool)
            except Exception as exc:
                print(f"[blast] www poll error (will retry): {exc}", flush=True)
                time.sleep(poll)
                continue

            if status == "WAITING":
                print(
                    f"[blast] www still waiting on NCBI "
                    f"({elapsed:.0f}s / {effective_timeout:.0f}s RID={rid})",
                    flush=True,
                )
                time.sleep(poll)
                continue
            if status == "FAILED":
                msg = f"[blast] www search FAILED (RID={rid})\n{body[:500]}\n"
                print(msg, flush=True)
                return BlastResult(returncode=1, stdout=b"", stderr=msg.encode())
            if status == "UNKNOWN":
                msg = f"[blast] www unknown status (RID={rid})\n{body[:500]}\n"
                print(msg, flush=True)
                return BlastResult(returncode=1, stdout=b"", stderr=msg.encode())

            try:
                n = _xml_to_blast_tsv(body, out)
            except Exception as exc:
                msg = f"[blast] www XML→TSV failed: {exc}\n"
                print(msg, flush=True)
                return BlastResult(returncode=1, stdout=b"", stderr=msg.encode())
            print(
                f"[blast] www finished RID={rid} hits_rows={n} out={out} "
                f"elapsed={time.time() - started:.0f}s",
                flush=True,
            )
            return BlastResult(returncode=0, stdout=b"", stderr=b"")
    finally:
        stop.set()
    # Defensive: try/finally without an explicit return yields None in some paths.
    return BlastResult(
        returncode=1,
        stdout=b"",
        stderr=b"[blast] www backend exited without a result\n",
    )


def _www_put(
    fasta: str,
    *,
    database: str,
    hitlist_size: int,
    word_size: int,
    evalue: float,
    email: str,
    tool: str,
) -> tuple[str, int]:
    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": database,
        "QUERY": fasta,
        "EXPECT": str(evalue),
        "WORD_SIZE": str(word_size),
        "HITLIST_SIZE": str(hitlist_size),
        "FILTER": "F",
        "FORMAT_TYPE": "XML",
        "CLIENT": "web",
        "EMAIL": email,
        "TOOL": tool,
    }
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(
        NCBI_BLAST_URL,
        data=data,
        method="POST",
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    rid_m = re.search(r"RID\s*=\s*([A-Za-z0-9_-]+)", text)
    rtoe_m = re.search(r"RTOE\s*=\s*(\d+)", text)
    if not rid_m:
        raise RuntimeError(f"no RID in PUT response: {text[:400]}")
    return rid_m.group(1), int(rtoe_m.group(1)) if rtoe_m else 30


def _www_get(rid: str, *, email: str, tool: str) -> tuple[str, str]:
    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "XML",
        "EMAIL": email,
        "TOOL": tool,
    }
    url = NCBI_BLAST_URL + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, method="GET")
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode("utf-8", errors="replace")

    if "Status=WAITING" in text or "Status = WAITING" in text:
        return "WAITING", text
    if "Status=FAILED" in text or "Status = FAILED" in text:
        return "FAILED", text
    if "Status=UNKNOWN" in text or "Status = UNKNOWN" in text:
        return "UNKNOWN", text
    # Ready results are BlastOutput XML (no Status= line, or Status=READY).
    if "<BlastOutput" in text or "Status=READY" in text or "Status = READY" in text:
        return "READY", text
    # Sometimes NCBI still returns HTML chrome while waiting.
    if "RID =" in text and "QBlastInfoBegin" in text:
        return "WAITING", text
    return "READY", text


def _xml_to_blast_tsv(xml_text: str, out_path: str) -> int:
    """Convert NCBI BlastOutput XML into PC ``BLAST_OUTFMT`` tabular TSV."""
    # Strip any leading HTML/status wrapper if present.
    start = xml_text.find("<BlastOutput")
    if start < 0:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        Path(out_path).write_text("")
        return 0
    xml_text = xml_text[start:]
    root = ET.fromstring(xml_text)

    query_id = ""
    for el in root.iter():
        if el.tag == "BlastOutput_query-def" or el.tag.endswith("}BlastOutput_query-def"):
            query_id = (el.text or "").split()[0] if el.text else ""
            break
        if el.tag == "Iteration_query-def" or el.tag.endswith("}Iteration_query-def"):
            query_id = (el.text or "").split()[0] if el.text else ""
            break

    rows: list[list[str]] = []
    for hit in root.iter():
        if not (hit.tag == "Hit" or hit.tag.endswith("}Hit")):
            continue
        hit_id = ""
        hit_acc = ""
        hit_def = ""
        for child in list(hit):
            tag = child.tag.split("}")[-1]
            if tag == "Hit_id":
                hit_id = (child.text or "").strip()
            elif tag == "Hit_accession":
                hit_acc = (child.text or "").strip()
            elif tag == "Hit_def":
                hit_def = (child.text or "").strip()
        if not hit_acc and hit_id:
            # ref|XP_123| or gi|…|ref|XP_123.1|
            m = re.search(r"(?:ref\|)?([A-Z0-9_]+)(?:\.\d+)?\|?", hit_id)
            if m:
                hit_acc = m.group(1)

        for hsp in hit.iter():
            if not (hsp.tag == "Hsp" or hsp.tag.endswith("}Hsp")):
                continue
            # Skip nested false positives: only real Hsp elements have Hsp_* kids.
            hsp_map = {}
            for c in list(hsp):
                hsp_map[c.tag.split("}")[-1]] = (c.text or "").strip()
            if "Hsp_evalue" not in hsp_map and "Hsp_bit-score" not in hsp_map:
                continue

            fields = {name: "" for name in BLAST_TSV_FIELDS}
            fields["qseqid"] = query_id
            fields["sseqid"] = hit_id
            fields["sacc"] = hit_acc
            fields["saccver"] = hit_acc
            fields["scomnames"] = hit_def[:80].replace("\t", " ")

            align_len = float(hsp_map.get("Hsp_align-len") or 0) or 1.0
            identity = float(hsp_map.get("Hsp_identity") or 0)
            fields["pident"] = f"{(100.0 * identity / align_len):.3f}"
            fields["length"] = hsp_map.get("Hsp_align-len", "")
            fields["mismatch"] = str(
                max(
                    0,
                    int(align_len)
                    - int(float(hsp_map.get("Hsp_identity") or 0))
                    - int(float(hsp_map.get("Hsp_gaps") or 0)),
                )
            )
            fields["gapopen"] = hsp_map.get("Hsp_gaps", "0")
            fields["qstart"] = hsp_map.get("Hsp_query-from", "")
            fields["qend"] = hsp_map.get("Hsp_query-to", "")
            fields["sstart"] = hsp_map.get("Hsp_hit-from", "")
            fields["send"] = hsp_map.get("Hsp_hit-to", "")
            fields["evalue"] = hsp_map.get("Hsp_evalue", "")
            fields["bitscore"] = hsp_map.get("Hsp_bit-score", "")
            rows.append([fields[name] for name in BLAST_TSV_FIELDS])

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for row in rows:
            fh.write("\t".join(row) + "\n")
    return len(rows)


# --------------------------------------------------------------------------- #
# Legacy blastp -remote (timeout-guarded)
# --------------------------------------------------------------------------- #


def _run_blast_remote(
    query: str,
    out: str,
    max_target_seqs: int,
    outfmt: str,
    word_size: int,
    evalue: float,
) -> BlastResult:
    database = os.environ.get("PC_BLAST_DB", "nr").strip() or "nr"
    timeout = float(os.environ.get("PC_BLAST_TIMEOUT", "900"))
    cmd = [
        "blastp",
        "-remote",
        "-db",
        database,
        "-query",
        query,
        "-out",
        out,
        "-max_target_seqs",
        str(max_target_seqs),
        "-outfmt",
        outfmt,
        "-word_size",
        str(word_size),
        "-evalue",
        str(evalue),
    ]
    print(
        f"[blast] remote blastp db={database} "
        f"max_target_seqs={max_target_seqs} word_size={word_size} "
        f"timeout={timeout:.0f}s query={query}",
        flush=True,
    )

    stop = threading.Event()

    def _heartbeat():
        started = time.time()
        while not stop.wait(30.0):
            elapsed = time.time() - started
            print(
                f"[blast] still waiting on NCBI remote blastp ({elapsed:.0f}s / {timeout:.0f}s)",
                flush=True,
            )

    hb = threading.Thread(target=_heartbeat, name="blast-remote-hb", daemon=True)
    hb.start()
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
        )
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            stop.set()
            try:
                os.killpg(proc.pid, 15)
            except ProcessLookupError:
                pass
            try:
                stdout, stderr = proc.communicate(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, 9)
                except ProcessLookupError:
                    pass
                stdout, stderr = proc.communicate()
            msg = (
                f"\n[blast] TIMEOUT after {timeout:.0f}s waiting on NCBI remote blastp\n"
            ).encode()
            print(f"[blast] TIMEOUT after {timeout:.0f}s", flush=True)
            return BlastResult(
                returncode=124,
                stdout=stdout or b"",
                stderr=(stderr or b"") + msg,
            )
        result = BlastResult(returncode=proc.returncode, stdout=stdout, stderr=stderr)
    finally:
        stop.set()

    if result.stderr:
        try:
            err_tail = result.stderr.decode("utf-8", errors="replace")[-500:]
        except Exception:
            err_tail = "<binary stderr>"
        if err_tail.strip():
            print(f"[blast] blastp stderr (tail):\n{err_tail}", flush=True)
    print(f"[blast] blastp finished rc={result.returncode}", flush=True)

    # A refused request exits zero with the error only on stderr, so the exit code is reported
    # back as a non-failure and the caller's word-size backoff never runs. Surface it as the
    # failure it is, keeping the original stderr for the caller's error message.
    if blast_call_failed(result) and result.returncode == 0:
        print("[blast] blastp reported an error despite exiting zero; treating as a failure")
        return BlastResult(returncode=1, stdout=result.stdout, stderr=result.stderr)

    return result
