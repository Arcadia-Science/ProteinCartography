#!/usr/bin/env python
"""Map RefSeq / EMBL CDS accessions to UniProtKB.

Handles UniProt idmapping status correctly (303 redirect / FINISHED / ERROR),
batches large ID lists, and retries transient ``jobStatus: ERROR`` responses.
"""

from __future__ import annotations
import argparse
import os
import sys
import time
from time import sleep

import pandas as pd
from api_utils import (
    UniProtWithExpBackoff,
    session_with_retry,
)
from constants import UniProtService

# If necessary, mock bioservices mapping (see ``tests.mocks``).
if os.environ.get("PROTEINCARTOGRAPHY_WAS_CALLED_BY_PYTEST") == "true":
    from tests import mocks

    mocks.mock_bioservices_uniprot_mapping()

__all__ = ["map_refseqids_bioservices", "map_refseqids_rest"]

DEFAULT_DBS = ["EMBL-GenBank-DDBJ_CDS", "RefSeq_Protein"]
UNIPROT_IDMAPPING_API = "https://rest.uniprot.org/idmapping"

# Per-batch polling. Override via env.
BATCH_SIZE = int(os.environ.get("PC_UNIPROT_MAP_BATCH", "500"))
MAX_WAIT_SECONDS = int(os.environ.get("PC_UNIPROT_MAP_TIMEOUT", "1800"))
POLL_SECONDS = float(os.environ.get("PC_UNIPROT_MAP_POLL", "10"))
SUBMIT_RETRIES = int(os.environ.get("PC_UNIPROT_MAP_RETRIES", "5"))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-d", "--databases", nargs="+", default=DEFAULT_DBS)
    parser.add_argument("-s", "--service", default=UniProtService.REST.value)
    return parser.parse_args()


def map_refseqids_bioservices(
    input_file: str, output_file: str, query_dbs: list, return_full=False
):
    uniprot = UniProtWithExpBackoff()
    with open(input_file) as f:
        ids = f.read().splitlines()
    ids = ids[:100000]

    dummy_df = pd.DataFrame()
    for i, db in enumerate(query_dbs):
        results = uniprot.mapping(
            db, "UniProtKB", query=",".join(ids), max_waiting_time=MAX_WAIT_SECONDS
        )
        if not results or "results" not in results:
            continue
        results_df = pd.json_normalize(results["results"])
        if len(results_df) == 0:
            continue
        dummy_df = (
            results_df if i == 0 or dummy_df.empty else pd.concat([dummy_df, results_df], axis=0)
        )

    hits = (
        dummy_df["to.primaryAccession"].unique()
        if not dummy_df.empty and "to.primaryAccession" in dummy_df.columns
        else []
    )
    with open(output_file, "w+") as f:
        f.writelines(hit + "\n" for hit in hits)
    if return_full:
        return dummy_df


def _chunks(items: list[str], size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _poll_job(session, job_id: str) -> dict:
    """Wait until UniProt finishes (303 / FINISHED / results) or errors."""
    deadline = time.time() + MAX_WAIT_SECONDS
    last_status = None
    while time.time() < deadline:
        response = session.get(
            f"{UNIPROT_IDMAPPING_API}/status/{job_id}",
            allow_redirects=False,
        )
        # Finished jobs redirect to the results collection.
        if response.status_code in (303, 302):
            return {"jobStatus": "FINISHED"}
        try:
            payload = response.json()
        except ValueError:
            payload = {}
        last_status = payload
        status = payload.get("jobStatus")
        if status == "FINISHED" or "results" in payload:
            return payload if "results" in payload else {"jobStatus": "FINISHED"}
        if status == "ERROR":
            return payload
        sleep(POLL_SECONDS)
    raise TimeoutError(
        f"UniProt idmapping job {job_id} did not finish within {MAX_WAIT_SECONDS}s; "
        f"last status={last_status}"
    )


def _fetch_results(session, job_id: str) -> list[dict]:
    stream = session.get(f"{UNIPROT_IDMAPPING_API}/stream/{job_id}")
    stream.raise_for_status()
    payload = stream.json()
    return list(payload.get("results") or [])


def _map_batch(session, db: str, batch: list[str]) -> list[dict]:
    """Submit one batch with retries on UniProt-side ERROR / transport flakes."""
    last_error = None
    for attempt in range(1, SUBMIT_RETRIES + 1):
        try:
            ticket = session.post(
                f"{UNIPROT_IDMAPPING_API}/run",
                {"ids": ",".join(batch), "from": db, "to": "UniProtKB"},
            )
            ticket.raise_for_status()
            job_id = ticket.json()["jobId"]
            status = _poll_job(session, job_id)
            if status.get("jobStatus") == "ERROR":
                last_error = status.get("errors") or status
                print(
                    f"[map_refseq_ids] UniProt ERROR on {db} batch "
                    f"(n={len(batch)}, attempt {attempt}/{SUBMIT_RETRIES}): {last_error}",
                    file=sys.stderr,
                )
                sleep(min(60, POLL_SECONDS * attempt * 2))
                continue
            if "results" in status:
                return list(status.get("results") or [])
            return _fetch_results(session, job_id)
        except Exception as exc:  # noqa: BLE001 — retry transient API/transport failures
            last_error = exc
            print(
                f"[map_refseq_ids] {type(exc).__name__} on {db} batch "
                f"(n={len(batch)}, attempt {attempt}/{SUBMIT_RETRIES}): {exc}",
                file=sys.stderr,
            )
            sleep(min(60, POLL_SECONDS * attempt * 2))
    raise RuntimeError(
        f"UniProt idmapping failed for {db} after {SUBMIT_RETRIES} attempts "
        f"(batch size {len(batch)}): {last_error}"
    )


def map_refseqids_rest(input_file: str, output_file: str, query_dbs: list, return_full=False):
    with open(input_file) as f:
        input_ids = list(dict.fromkeys(line.strip() for line in f if line.strip()))

    if not input_ids:
        print("[map_refseq_ids] empty input — writing empty UniProt hit list", file=sys.stderr)
        with open(output_file, "w+") as f:
            pass
        if return_full:
            return pd.DataFrame()
        return

    session = session_with_retry()
    frames: list[pd.DataFrame] = []

    for db in query_dbs:
        print(
            f"[map_refseq_ids] mapping {len(input_ids)} ids via {db} "
            f"(batch={BATCH_SIZE}, timeout={MAX_WAIT_SECONDS}s)",
            file=sys.stderr,
        )
        rows: list[dict] = []
        for batch in _chunks(input_ids, max(1, BATCH_SIZE)):
            rows.extend(_map_batch(session, db, batch))
        if not rows:
            continue
        frames.append(pd.DataFrame(rows))

    dummy_df = pd.concat(frames, axis=0) if frames else pd.DataFrame()
    if dummy_df.empty:
        hits = []
    elif "to" in dummy_df.columns:
        hits = dummy_df["to"].dropna().unique()
    elif "to.primaryAccession" in dummy_df.columns:
        hits = dummy_df["to.primaryAccession"].dropna().unique()
    else:
        hits = []

    with open(output_file, "w+") as f:
        f.writelines(f"{hit}\n" for hit in hits)

    if return_full:
        return dummy_df


def main():
    args = parse_args()
    service = UniProtService(args.service)
    if service == UniProtService.BIOSERVICES:
        map_refseqids_bioservices(args.input, args.output, args.databases)
    elif service == UniProtService.REST:
        map_refseqids_rest(args.input, args.output, args.databases)
    else:
        sys.exit(f"unknown UniProt service: {args.service}")


if __name__ == "__main__":
    main()
