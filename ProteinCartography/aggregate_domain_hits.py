#!/usr/bin/env python
"""Union domain-path BLAST/Foldseek hit lists and record which query domain found each hit."""

from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path

import domain_utils as du

__all__ = ["aggregate_domain_hits"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Hit files named {domain_id}.blast_hits.uniprot.txt or {domain_id}.foldseek_hits.txt",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Aggregated unique accessions, one per line.",
    )
    parser.add_argument(
        "-f",
        "--found-by",
        required=True,
        help="TSV of accession and found_by (comma-separated query domain ids).",
    )
    return parser.parse_args()


def _domain_id_from_hit_file(path: str) -> str:
    stem = Path(path).name
    for suffix in (
        ".blast_hits.uniprot.txt",
        ".foldseek_hits.txt",
        ".blast_hits.refseq.txt",
    ):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return Path(path).stem


def aggregate_domain_hits(input_files: list, output_file: str, found_by_file: str) -> list[str]:
    found_by: dict[str, set[str]] = defaultdict(set)
    for file in input_files:
        domain_id = _domain_id_from_hit_file(file)
        try:
            du.parse_domain_id(domain_id)
        except ValueError:
            domain_id = Path(file).stem
        with open(file) as handle:
            for line in handle:
                acc = line.strip()
                if not acc:
                    continue
                found_by[acc].add(domain_id)

    accessions = sorted(found_by)
    Path(output_file).write_text("".join(acc + "\n" for acc in accessions))
    found_by_path = Path(found_by_file)
    found_by_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["accession\tfound_by\n"]
    for acc in accessions:
        sources = ",".join(sorted(found_by[acc]))
        lines.append(f"{acc}\t{sources}\n")
    found_by_path.write_text("".join(lines))
    return accessions


def main():
    args = parse_args()
    aggregate_domain_hits(args.input, args.output, args.found_by)


if __name__ == "__main__":
    main()
