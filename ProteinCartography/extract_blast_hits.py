#!/usr/bin/env python
import argparse

import constants
import pandas as pd

# only import these functions when using import *
__all__ = ["extract_blast_hits"]


# parse command line arguments
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
    args = parser.parse_args()

    return args


# take an input blastresults file and create a .txt file from that
def extract_blast_hits(input_file: str, output_file: str, column_names: list):
    """
    Takes an input blast_results.tsv file, reads the accessions,
    and prints unique hits to a .txt file, one per line.

    Args:
        input_file (str): path of input blast_results.tsv file.
        output_file (str): path of destination blast_hits.txt file.
        names (str): names of columns of blast results.
    """
    df = pd.read_csv(input_file, sep="\t", names=column_names)

    hits = df["sacc"].unique()

    if len(hits) == 0:
        raise Exception("No hits were returned. Check to see if remote BLAST failed.")

    with open(output_file, "w+") as f:
        f.writelines(hit + "\n" for hit in hits)


# run this if called from the interpreter
def main():
    args = parse_args()
    blast_column_names = [name for name in args.blast_format_string.split(" ") if name != "6"]
    extract_blast_hits(
        input_file=args.input, output_file=args.output, column_names=blast_column_names
    )


# check if called from interpreter
if __name__ == "__main__":
    main()
