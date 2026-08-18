#!/usr/bin/env python
"""Assign TED/user domain boundaries and crop query or hit structures.

Query-gate: TED (or a user TSV) on query proteins only. The domain DAG runs
only when at least one query has two or more kept domains. TED failures on a
query never fail the protein pipeline.

Hit assignment: TED on domain-path hit accessions; misses are omitted.
"""

from __future__ import annotations
import argparse
import csv
import os
import sys
from pathlib import Path

import domain_utils as du

# If necessary, mock the web API responses (see ``tests.mocks``).
if os.environ.get("PROTEINCARTOGRAPHY_SHOULD_USE_MOCKS") == "true":
    from tests import mocks

    mocks.mock_requests_session_request()

__all__ = [
    "fetch_ted_summary",
    "load_user_domains",
    "assign_parent",
    "run_query_gate",
    "assign_and_crop_hits",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    gate = sub.add_parser("query-gate", help="Assign query domains and write the gate files.")
    gate.add_argument("--protids", nargs="+", required=True, help="Query protein ids.")
    gate.add_argument("--input-dir", required=True, help="Directory with query FASTA/PDB files.")
    gate.add_argument("--output-dir", required=True, help="Domain-path directory to write.")
    gate.add_argument("--user-domains-file", default="", help="Optional user domain TSV.")
    gate.add_argument(
        "--min-domain-length",
        type=int,
        default=du.DEFAULT_MIN_DOMAIN_LENGTH,
    )
    gate.add_argument(
        "--domain-map",
        default="auto",
        choices=["auto", "off"],
    )

    hits = sub.add_parser("assign-hits", help="TED-assign and crop domain-path hit PDBs.")
    hits.add_argument("--accessions-file", required=True, help="One UniProt accession per line.")
    hits.add_argument("--pdb-dir", required=True, help="Directory of downloaded full-length PDBs.")
    hits.add_argument("--output-dir", required=True, help="Directory for cropped domain PDBs.")
    hits.add_argument("--output-tsv", required=True, help="Path to domain_features.tsv.")
    hits.add_argument("--cache-dir", required=True, help="TED JSON cache directory.")
    hits.add_argument(
        "--min-domain-length",
        type=int,
        default=du.DEFAULT_MIN_DOMAIN_LENGTH,
    )
    hits.add_argument(
        "--query-domains-tsv",
        default="",
        help="Query domain rows to prepend so query domains are always on the map.",
    )
    hits.add_argument(
        "--query-structures-dir",
        default="",
        help="Directory of already-cropped query domain PDBs to copy into output-dir.",
    )
    hits.add_argument("--user-domains-file", default="")

    return parser.parse_args()


def _ted_rate_limited_get(session, url: str, timeout: int = 60):
    """GET with the same polite rate limit pattern as download_pdbs."""
    try:
        from ratelimiter import RateLimiter
    except ImportError:
        return session.get(url, timeout=timeout)

    limiter = getattr(_ted_rate_limited_get, "_limiter", None)
    if limiter is None:
        limiter = RateLimiter(max_calls=10, period=1)
        _ted_rate_limited_get._limiter = limiter
    with limiter:
        return session.get(url, timeout=timeout)


def _ted_get_page(session, acc: str, skip: int, limit: int) -> dict | None:
    """Fetch one TED summary page. 429 is retried; other errors return None."""
    url = du.TED_SUMMARY_URL.format(acc=acc) + f"?skip={skip}&limit={limit}"
    last_status = None
    for attempt in range(1, du.TED_MAX_429_RETRIES + 1):
        try:
            response = _ted_rate_limited_get(session, url, timeout=60)
        except Exception as exc:
            print(f"[assign_domains] TED request failed for {acc}: {exc}", file=sys.stderr)
            return None

        last_status = getattr(response, "status_code", None)
        if last_status == 429:
            print(
                f"[assign_domains] TED HTTP 429 for {acc} (attempt {attempt}/"
                f"{du.TED_MAX_429_RETRIES}); retrying",
                file=sys.stderr,
            )
            continue
        if last_status == 404:
            return {"data": [], "count": 0}
        if not getattr(response, "ok", False) or (last_status is not None and last_status >= 400):
            print(
                f"[assign_domains] TED HTTP {last_status} for {acc}; treating as miss",
                file=sys.stderr,
            )
            return None
        try:
            payload = response.json()
        except Exception as exc:
            print(f"[assign_domains] TED JSON error for {acc}: {exc}", file=sys.stderr)
            return None
        if not isinstance(payload, dict):
            payload = {"data": payload or [], "count": 0}
        return payload

    print(
        f"[assign_domains] TED HTTP {last_status} for {acc} after retries; treating as miss",
        file=sys.stderr,
    )
    return None


def fetch_ted_summary(accession: str, session, cache_dir: Path) -> dict | None:
    """Return a TED summary payload or None on miss/error. Caches successful GETs and 404s.

    Pages with skip/limit until a short page is returned. ``count`` on the API is the
    page size, not a total, so pagination stops on a short page rather than count.
    """
    acc = du.canonical_uniprot_accession(accession)
    cached = du.load_ted_cache(cache_dir, acc)
    if cached is not None:
        return cached

    all_data: list = []
    skip = 0
    while True:
        payload = _ted_get_page(session, acc, skip, du.TED_PAGE_SIZE)
        if payload is None:
            if skip == 0:
                return None
            break
        data = payload.get("data") or []
        all_data.extend(data)
        if len(data) < du.TED_PAGE_SIZE:
            break
        skip += du.TED_PAGE_SIZE

    merged = {"data": all_data, "count": len(all_data)}
    du.save_ted_cache(cache_dir, acc, merged)
    return merged


def load_user_domains(path: str) -> dict[str, list[dict]]:
    """Map parent_protid → domain rows from an optional user TSV."""
    if not path:
        return {}
    user_path = Path(path)
    if not user_path.is_file():
        return {}
    with user_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        records = list(reader)
    rows = du.rows_from_user_tsv_records(records)
    by_parent: dict[str, list[dict]] = {}
    for row in rows:
        by_parent.setdefault(row["parent_protid"], []).append(row)
    return by_parent


def assign_parent(
    parent_protid: str,
    session,
    cache_dir: Path,
    user_by_parent: dict[str, list[dict]],
    min_domain_length: int,
) -> list[dict]:
    """User TSV wins; else TED; else empty. Never raises."""
    if parent_protid in user_by_parent:
        return du.filter_domain_rows(user_by_parent[parent_protid], min_domain_length)

    if not du.looks_like_uniprot_accession(parent_protid):
        return []

    payload = fetch_ted_summary(parent_protid, session, cache_dir)
    if not payload:
        return []
    rows = du.rows_from_ted_payload(parent_protid, payload)
    return du.filter_domain_rows(rows, min_domain_length)


def _write_domains_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=du.DOMAIN_FEATURES_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in du.DOMAIN_FEATURES_COLUMNS})


def _crop_query_files(parent: str, rows: list[dict], input_dir: Path, query_structures: Path):
    fasta = input_dir / f"{parent}.fasta"
    pdb = input_dir / f"{parent}.pdb"
    for row in rows:
        domain_id = row["protid"]
        chopping = row["chopping"]
        if fasta.is_file():
            du.crop_fasta_file(fasta, chopping, query_structures / f"{domain_id}.fasta", domain_id)
        if pdb.is_file():
            dest = query_structures / f"{domain_id}.pdb"
            n_ca = du.crop_pdb_file(pdb, chopping, dest)
            row["nres_domain"] = n_ca
            if n_ca <= 0:
                dest.unlink(missing_ok=True)
                fasta_out = query_structures / f"{domain_id}.fasta"
                if fasta_out.is_file():
                    fasta_out.unlink()


def run_query_gate(
    protids: list[str],
    input_dir: str,
    output_dir: str,
    user_domains_file: str = "",
    min_domain_length: int = du.DEFAULT_MIN_DOMAIN_LENGTH,
    domain_map: str = "auto",
    session=None,
) -> str:
    """Write gate.txt, query_domains.tsv, query_domain_ids.txt, and cropped query files.

    Returns ``\"on\"`` or ``\"off\"``.
    """
    output = Path(output_dir)
    cache_dir = output / "ted_cache"
    query_structures = output / "query_structures"
    query_structures.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    gate_path = output / "gate.txt"
    ids_path = output / "query_domain_ids.txt"
    tsv_path = output / "query_domains.tsv"

    if str(domain_map).lower() == "off":
        gate_path.write_text("off\n")
        ids_path.write_text("")
        _write_domains_tsv(tsv_path, [])
        return "off"

    if session is None:
        from api_utils import session_with_retry

        session = session_with_retry()
    user_by_parent = load_user_domains(user_domains_file)

    multi_rows: list[dict] = []
    for parent in protids:
        rows = assign_parent(parent, session, cache_dir, user_by_parent, min_domain_length)
        if not du.is_multidomain(len(rows)):
            continue
        _crop_query_files(parent, rows, Path(input_dir), query_structures)
        cropped = [row for row in rows if int(row.get("nres_domain") or 0) > 0]
        if du.is_multidomain(len(cropped)):
            multi_rows.extend(cropped)

    gate = "on" if multi_rows else "off"
    gate_path.write_text(gate + "\n")
    ids_path.write_text("".join(row["protid"] + "\n" for row in multi_rows))
    _write_domains_tsv(tsv_path, multi_rows)
    return gate


def assign_and_crop_hits(
    accessions_file: str,
    pdb_dir: str,
    output_dir: str,
    output_tsv: str,
    cache_dir: str,
    min_domain_length: int = du.DEFAULT_MIN_DOMAIN_LENGTH,
    query_domains_tsv: str = "",
    query_structures_dir: str = "",
    user_domains_file: str = "",
    session=None,
) -> list[dict]:
    """TED-assign hit PDBs, omit misses, crop survivors, prepend query domain rows."""
    if session is None:
        from api_utils import session_with_retry

        session = session_with_retry()
    user_by_parent = load_user_domains(user_domains_file)
    pdb_dir_path = Path(pdb_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    accessions = [
        line.strip() for line in Path(accessions_file).read_text().splitlines() if line.strip()
    ]

    query_rows: list[dict] = []
    query_parents: set[str] = set()
    if query_domains_tsv and Path(query_domains_tsv).is_file():
        with Path(query_domains_tsv).open(newline="") as handle:
            query_rows = list(csv.DictReader(handle, delimiter="\t"))
        query_parents = {
            str(row.get("parent_protid", "")).strip()
            for row in query_rows
            if str(row.get("parent_protid", "")).strip()
        }
        query_src = Path(query_structures_dir) if query_structures_dir else None
        for row in query_rows:
            dest = out_dir / f"{row['protid']}.pdb"
            src = (query_src / f"{row['protid']}.pdb") if query_src else None
            if src is not None and src.is_file() and not dest.exists():
                dest.write_text(src.read_text())

    hit_rows: list[dict] = []
    for acc in accessions:
        if acc in query_parents:
            continue
        rows = assign_parent(acc, session, cache, user_by_parent, min_domain_length)
        if not rows:
            print(f"[assign_domains] omitting hit {acc}: no TED/user domains", file=sys.stderr)
            continue
        pdb_path = pdb_dir_path / f"{acc}.pdb"
        if not pdb_path.is_file():
            print(f"[assign_domains] omitting hit {acc}: missing PDB", file=sys.stderr)
            continue
        for row in rows:
            dest = out_dir / f"{row['protid']}.pdb"
            n_ca = du.crop_pdb_file(pdb_path, row["chopping"], dest)
            if n_ca <= 0:
                dest.unlink(missing_ok=True)
                print(
                    f"[assign_domains] omitting {row['protid']}: empty crop for {acc}",
                    file=sys.stderr,
                )
                continue
            row["nres_domain"] = n_ca
            hit_rows.append(row)

    all_rows = query_rows + hit_rows
    # DictReader values are strings; normalize types for the TSV writer.
    normalized = []
    for row in all_rows:
        normalized.append(
            {
                col: row.get(col, "")
                if col not in {"domain_index", "start", "end", "nres_domain"}
                else int(row[col])
                if str(row.get(col, "")).strip()
                else ""
                for col in du.DOMAIN_FEATURES_COLUMNS
            }
        )
    _write_domains_tsv(Path(output_tsv), normalized)
    return normalized


def main():
    args = parse_args()
    if args.command == "query-gate":
        run_query_gate(
            protids=args.protids,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            user_domains_file=args.user_domains_file,
            min_domain_length=args.min_domain_length,
            domain_map=args.domain_map,
        )
    elif args.command == "assign-hits":
        assign_and_crop_hits(
            accessions_file=args.accessions_file,
            pdb_dir=args.pdb_dir,
            output_dir=args.output_dir,
            output_tsv=args.output_tsv,
            cache_dir=args.cache_dir,
            min_domain_length=args.min_domain_length,
            query_domains_tsv=args.query_domains_tsv,
            query_structures_dir=args.query_structures_dir,
            user_domains_file=args.user_domains_file,
        )


if __name__ == "__main__":
    main()
