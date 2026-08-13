#!/usr/bin/env python
"""Fetch UniProt metadata for aggregated ProteinCartography hits.

Uses UniProt ``/uniprotkb/accessions`` (REST) to avoid the 100-OR search limit,
with retries and safe handling of empty/null responses.
"""

from __future__ import annotations
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from api_utils import UniProtWithExpBackoff, session_with_retry
from constants import UniProtService

# If necessary, mock bioservices search (see ``tests.mocks``).
if os.environ.get("PROTEINCARTOGRAPHY_SHOULD_USE_MOCKS") == "true":
    from tests import mocks

    mocks.mock_bioservices_uniprot_search()

__all__ = ["query_uniprot"]

REQUIRED_FIELDS_DICT = {
    "Entry": "accession",
    "Entry Name": "id",
    "Protein names": "protein_name",
    "Gene Names (primary)": "gene_primary",
    "Annotation": "annotation_score",
    "Organism": "organism_name",
    "Taxonomic lineage": "lineage",
    "Length": "length",
    "Fragment": "fragment",
    "Sequence": "sequence",
}
OTHER_FIELDS_DICT = {
    "Reviewed": "reviewed",
    "Gene Names": "gene_names",
    "Protein existence": "protein_existence",
    "Sequence version": "sequence_version",
    "RefSeq": "xref_refseq",
    "GeneID": "xref_geneid",
    "EMBL": "xref_embl",
    "AlphaFoldDB": "xref_alphafolddb",
    "PDB": "xref_pdb",
    "Pfam": "xref_pfam",
    "InterPro": "xref_interpro",
}
DEFAULT_FIELDS_DICT = REQUIRED_FIELDS_DICT | OTHER_FIELDS_DICT
REQUIRED_FIELDS = list(REQUIRED_FIELDS_DICT.values())
DEFAULT_FIELDS = list(DEFAULT_FIELDS_DICT.values())

UNIPROT_ACCESSIONS = "https://rest.uniprot.org/uniprotkb/accessions"
BATCH_SIZE = int(os.environ.get("PC_UNIPROT_META_BATCH", "100"))
RETRIES = int(os.environ.get("PC_UNIPROT_META_RETRIES", "5"))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument(
        "-s",
        "--service",
        default=UniProtService.REST.value,
        help="How to fetch metadata (rest|bioservices). REST uses /uniprotkb/accessions.",
    )
    parser.add_argument(
        "-a",
        "--additional-fields",
        nargs="*",
        help="additional non-default fields to fetch from uniprot if using REST",
    )
    return parser.parse_args()


def _chunks(items: list[str], size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _fetch_accessions_tsv(session, accessions: list[str], fields_string: str) -> str:
    """GET /uniprotkb/accessions — avoids the 100-OR search limit."""
    last_error = None
    for attempt in range(1, RETRIES + 1):
        try:
            response = session.get(
                UNIPROT_ACCESSIONS,
                params={
                    "accessions": ",".join(accessions),
                    "format": "tsv",
                    "fields": fields_string,
                },
                timeout=120,
            )
            if response.status_code == 400:
                # Invalid ids in the batch — shrink and retry halves if needed.
                last_error = f"400 Bad Request: {response.text[:200]}"
                print(
                    f"[fetch_uniprot_metadata] Bad Request on {len(accessions)} ids "
                    f"(attempt {attempt}/{RETRIES}): {last_error}",
                    file=sys.stderr,
                )
                if len(accessions) > 1 and attempt < RETRIES:
                    mid = len(accessions) // 2
                    left = _fetch_accessions_tsv(session, accessions[:mid], fields_string)
                    right = _fetch_accessions_tsv(session, accessions[mid:], fields_string)
                    # Merge TSVs (drop duplicate header from right).
                    left_lines = left.splitlines()
                    right_lines = right.splitlines()
                    if not left_lines:
                        return right
                    if not right_lines:
                        return left
                    return "\n".join(left_lines + right_lines[1:]) + "\n"
                time.sleep(min(60, 2 * attempt))
                continue
            response.raise_for_status()
            if not response.text or not response.text.strip():
                raise RuntimeError("empty UniProt accessions response")
            return response.text
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            print(
                f"[fetch_uniprot_metadata] {type(exc).__name__} "
                f"(n={len(accessions)}, attempt {attempt}/{RETRIES}): {exc}",
                file=sys.stderr,
            )
            time.sleep(min(60, 2 * attempt))
    raise RuntimeError(
        f"UniProt accessions fetch failed after {RETRIES} attempts "
        f"(batch size {len(accessions)}): {last_error}"
    )


def query_uniprot(
    input_file: str,
    output_file: str,
    batch_size=None,
    sub_batch_size=100,
    fmt="tsv",
    fields=DEFAULT_FIELDS,
    service=UniProtService.REST,
):
    """Download UniProt metadata for one accession per input line."""
    batch_size = int(batch_size or BATCH_SIZE)
    temp_filepath = output_file + ".temp"

    if os.path.exists(output_file):
        print("Output file already exists at this location. Aborting.")
        return

    header_written = False
    existing_data = set()
    if os.path.exists(temp_filepath):
        try:
            existing_df = pd.read_csv(temp_filepath, sep="\t", usecols=["Entry"])
            if "Entry" in existing_df.columns:
                existing_data = set(existing_df["Entry"].values)
                print(f"Loaded {len(existing_data)} entries from {temp_filepath}")
                header_written = True
            else:
                print(f"No 'Entry' column found in {temp_filepath}. Skipping existing data.")
        except pd.errors.EmptyDataError:
            print(f"{temp_filepath} is empty. Skipping existing data loading.")

    fields_string = ",".join(fields)
    session = session_with_retry()

    with open(input_file) as q:
        query_accessions = [line.rstrip("\n").strip() for line in q.readlines() if line.strip()]
    query_accessions = [e for e in query_accessions if e not in existing_data]
    # De-dupe while preserving order.
    query_accessions = list(dict.fromkeys(query_accessions))

    accession_batches = list(_chunks(query_accessions, max(1, batch_size)))
    total = len(query_accessions)
    print(
        f"[fetch_uniprot_metadata] fetching {total} accessions via {service} (batch={batch_size})",
        file=sys.stderr,
    )

    with open(temp_filepath, "a") as temp_file:
        for i, accession_batch in enumerate(accession_batches):
            print(f">> Starting batch {i + 1} of {len(accession_batches)}")
            if service == UniProtService.REST:
                batch_text = _fetch_accessions_tsv(session, accession_batch, fields_string)
            elif service == UniProtService.BIOSERVICES:
                # Keep bioservices available but cap OR count at UniProt's limit.
                uniprot = UniProtWithExpBackoff()
                or_batches = list(_chunks(accession_batch, 100))
                parts = []
                for sub in or_batches:
                    query_string = "(" + " OR ".join(f"accession:{a}" for a in sub) + ")"
                    results = uniprot.search(
                        query_string, columns=fields_string, size=sub_batch_size, progress=False
                    )
                    if results is None:
                        raise RuntimeError(
                            f"bioservices UniProt.search returned None for {len(sub)} ids"
                        )
                    parts.append(results)
                # Merge TSV parts
                lines_out = []
                for idx, part in enumerate(parts):
                    plines = part.splitlines()
                    if not plines:
                        continue
                    if idx == 0:
                        lines_out.extend(plines)
                    else:
                        lines_out.extend(plines[1:])
                batch_text = "\n".join(lines_out) + ("\n" if lines_out else "")
            else:
                raise ValueError(f"Unknown service {service}")

            if not batch_text or not batch_text.strip():
                print(f"empty batch {i + 1}; skipping", file=sys.stderr)
                continue

            lines = batch_text.splitlines()
            if not header_written:
                print_lines = lines
                header_written = True
            else:
                print_lines = lines[1:]

            for line in print_lines:
                print(line, file=temp_file)

            progress = max(0, len(lines) - 1)
            print(f"downloaded {progress} hits for batch {i + 1} (of {total} total)")

    if not os.path.exists(temp_filepath) or os.path.getsize(temp_filepath) == 0:
        raise RuntimeError(f"no UniProt metadata written for {input_file}")

    df = pd.read_csv(temp_filepath, sep="\t")
    if "Entry" not in df.columns:
        raise RuntimeError(f"UniProt TSV missing Entry column; columns={list(df.columns)}")
    df.insert(0, "protid", df["Entry"].values)

    def lineage_string_splitter(lineage_string):
        if lineage_string is np.nan:
            return np.nan
        return [rank.split(" (")[0] for rank in str(lineage_string).split(", ")]

    try:
        df["Lineage"] = df["Taxonomic lineage"].apply(lineage_string_splitter)
    except Exception:
        pass

    df.to_csv(output_file, sep="\t", index=None)
    os.remove(temp_filepath)
    return df


def main():
    args = parse_args()
    service = UniProtService(args.service)
    fields = list(DEFAULT_FIELDS)
    if args.additional_fields is not None:
        fields += args.additional_fields
    query_uniprot(args.input, args.output, fields=fields, service=service)


if __name__ == "__main__":
    main()
