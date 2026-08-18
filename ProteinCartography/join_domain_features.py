#!/usr/bin/env python
"""Join UniProt protein metadata onto domain rows via parent_protid.

Domain ``protid`` stays the domain id. UniProt columns (Protein names, Organism,
...) are copied from the parent accession so existing plot/semantic scripts work.
"""

from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd

__all__ = ["join_domain_uniprot"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain-features", required=True)
    parser.add_argument("--uniprot-features", required=True)
    parser.add_argument("--found-by", default="", help="Optional accession\\tfound_by TSV.")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def join_domain_uniprot(
    domain_features: str,
    uniprot_features: str,
    output_file: str,
    found_by_file: str = "",
) -> pd.DataFrame:
    domain_df = pd.read_csv(domain_features, sep="\t")
    uniprot_path = Path(uniprot_features)
    if uniprot_path.is_file() and uniprot_path.stat().st_size > 0:
        uniprot_df = pd.read_csv(uniprot_path, sep="\t")
        if "protid" in uniprot_df.columns:
            uniprot_df = uniprot_df.rename(columns={"protid": "parent_protid"})
            overlap = [
                c for c in uniprot_df.columns if c in domain_df.columns and c != "parent_protid"
            ]
            uniprot_df = uniprot_df.drop(columns=overlap)
            domain_df = domain_df.merge(uniprot_df, on="parent_protid", how="left")

    if found_by_file and Path(found_by_file).is_file():
        found_df = pd.read_csv(found_by_file, sep="\t")
        if "accession" in found_df.columns:
            found_df = found_df.rename(columns={"accession": "parent_protid"})
            domain_df = domain_df.merge(found_df, on="parent_protid", how="left")

    # Domain map Length is the cropped domain, not the parent chain.
    if "nres_domain" in domain_df.columns:
        domain_df["Length"] = domain_df["nres_domain"]

    domain_df.to_csv(output_file, sep="\t", index=False)
    return domain_df


def main():
    args = parse_args()
    join_domain_uniprot(
        args.domain_features, args.uniprot_features, args.output, found_by_file=args.found_by
    )


if __name__ == "__main__":
    main()
