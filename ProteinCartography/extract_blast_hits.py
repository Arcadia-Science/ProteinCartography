#!/usr/bin/env python
"""Extract RefSeq accessions from a BLAST results TSV.

An empty BLAST results file writes an empty hits file instead of raising, so
Foldseek-only maps can proceed when BLAST soft-fails or returns no hits.
"""

from __future__ import annotations
import argparse
import os
from pathlib import Path

import constants
import pandas as pd

__all__ = ["extract_blast_hits"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", required=True, help="path of input blast_results.tsv file."
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help=(
            "path of destination blast_hits.txt file, where blast hit accessions will be printed, "
            "one per line."
        ),
    )
    parser.add_argument(
        "-B",
        "--blast-format-string",
        default=constants.BLAST_OUTFMT,
        help=f"BLAST query format string.\n Defaults to '{constants.BLAST_OUTFMT}'",
    )
    return parser.parse_args()


def extract_blast_hits(input_file: str, output_file: str, column_names: list):
    """
    Takes an input blast_results.tsv file, reads the accessions,
    and prints unique hits to a .txt file, one per line.
    """
    path = Path(input_file)
    # Empty hits are allowed by default. Set PC_BLAST_SOFT_FAIL=0 to restore the
    # historical hard failure on empty BLAST results.
    soft = os.environ.get("PC_BLAST_SOFT_FAIL", "1").strip().lower() not in {
        "0",
        "false",
        "no",
        "off",
    }

    if not path.exists() or path.stat().st_size == 0:
        print(f"[blast] no BLAST results in {input_file}; writing empty hits file", flush=True)
        Path(output_file).write_text("")
        if soft:
            return
        raise Exception("No hits were returned. Check to see if remote BLAST failed.")

    df = pd.read_csv(input_file, sep="\t", names=column_names)
    if "sacc" not in df.columns or df.empty:
        print("[blast] BLAST results missing sacc / empty; writing empty hits file", flush=True)
        Path(output_file).write_text("")
        if soft:
            return
        raise Exception("No hits were returned. Check to see if remote BLAST failed.")

    hits = df["sacc"].dropna().unique()
    if len(hits) == 0:
        print("[blast] zero sacc hits; writing empty hits file", flush=True)
        Path(output_file).write_text("")
        if soft:
            return
        raise Exception("No hits were returned. Check to see if remote BLAST failed.")

    with open(output_file, "w+") as f:
        f.writelines(hit + "\n" for hit in hits)


def main():
    args = parse_args()
    blast_column_names = [name for name in args.blast_format_string.split(" ") if name != "6"]
    extract_blast_hits(
        input_file=args.input, output_file=args.output, column_names=blast_column_names
    )


if __name__ == "__main__":
    main()
