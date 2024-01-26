#!/usr/bin/env python
import argparse

import pandas as pd

__all__ = ["filter_results"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", required=True, help="path of input uniprot_features.tsv file"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="path of destination filtered hits list file",
    )
    parser.add_argument(
        "-m",
        "--min-length",
        default="0",
        help="minimum protein length of proteins to keep",
    )
    parser.add_argument(
        "-M",
        "--max-length",
        default="0",
        help="maximum protein length of proteins to keep. If set to 0, no upper limit is applied.",
    )
    parser.add_argument(
        "--excluded-protids",
        nargs="*",
        help="a list of protids to exclude from the results",
    )
    args = parser.parse_args()

    return args


def filter_results(
    input_file: str,
    output_file: str,
    filter_inactive=True,
    filter_fragment=True,
    min_length=0,
    max_length=0,
    excluded_protids=None,
):
    """
    Takes an input uniprot_features.tsv file and filters the results based on fragment status,
    inactive status, and size.

    Args:
        input_file (str): path of input tsv file with all uniprot features
        output_file (str): path of destination list text file where each accession is on a new line
        filter_inactive (bool): whether to remove inactive proteins
        filter_fragment (bool): whether to remove fragmentary proteins
        min_length (int): minimum length of proteins to keep
        max_length (int): maximum length of proteins to keep

    Returns:
        a pandas.DataFrame of the resulting features
    """

    df = pd.read_csv(input_file, sep="\t")

    filtered_df = df.copy(deep=True)

    if filter_inactive:
        filtered_df = filtered_df[
            (filtered_df["Protein names"] != "deleted") & ~filtered_df["Annotation"].isna()
        ]
    if filter_fragment:
        filtered_df = filtered_df[filtered_df["Fragment"] != "fragment"]

    if min_length > max_length and max_length != 0:
        raise ValueError("Minimum length must be less than maximum length.")

    if min_length > 0:
        filtered_df = filtered_df[filtered_df["Length"].astype(int) > min_length]
    if max_length > 0:
        filtered_df = filtered_df[filtered_df["Length"].astype(int) < max_length]

    if excluded_protids is None:
        excluded_protids = []

    filtered_df = filtered_df[~filtered_df["protid"].isin(excluded_protids)]

    with open(output_file, "w+") as f:
        f.writelines([protid + "\n" for protid in filtered_df["protid"]])

    return filtered_df


def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    excluded_protids = args.excluded_protids

    try:
        min_length = int(args.min_length)
    except (TypeError, ValueError):
        min_length = 0
    try:
        max_length = int(args.max_length)
    except (TypeError, ValueError):
        max_length = 0

    filter_results(
        input_file,
        output_file,
        min_length=min_length,
        max_length=max_length,
        excluded_protids=excluded_protids,
    )


if __name__ == "__main__":
    main()
