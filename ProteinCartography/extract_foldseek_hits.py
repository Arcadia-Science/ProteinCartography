#!/usr/bin/env python
import argparse
import os
import re

import constants
import pandas as pd

# only import these functions when using import *
__all__ = ["extract_foldseekhits"]

DEFAULT_EVALUE = 0.01


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", nargs="+", required=True, help="Takes .m8 file paths as input."
    )
    parser.add_argument("-o", "--output", required=True, help="Returns a .txt file as output.")
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=DEFAULT_EVALUE,
        help="Sets maximum evalue for filtering.",
    )
    parser.add_argument(
        "-m",
        "--max-num-hits",
        type=int,
        help=(
            "Maximum number of hits to include in the output file. "
            "If not provided, all hits will be included."
        ),
    )

    args = parser.parse_args()

    return args


def extract_foldseekhits(
    input_files: list, output_file: str, evalue=DEFAULT_EVALUE, max_num_hits=None
):
    """
    Takes a list of input tabular Foldseek results files from the API query (ending in .m8).
    Generates a .txt file containing a list of unique accessions across all the .m8 files.

    Args:
        input_files (list): list of string paths to input files.
        output_file (str): path of destination file.
    """

    # empty df for collecting results
    dummy_df = pd.DataFrame()

    # iterate through results files, reading them
    for i, file in enumerate(input_files):
        # load the file
        file_df = pd.read_csv(file, sep="\t", names=constants.FOLDSEEK_COLUMN_NAMES)

        if os.path.getsize(file) == 0:
            continue

        # extract the model ID from the results target column
        file_df["modelid"] = file_df["target"].str.split(" ", expand=True)[0]

        # extract only models that contain AF model string
        # this will need to be changed in the future
        file_df = file_df[file_df["modelid"].str.contains("-F1-model")]

        # filter by evalue
        file_df = file_df[file_df["evalue"] < evalue]

        # get the uniprot ID out from that target
        file_df["uniprotid"] = file_df["modelid"].apply(
            lambda x: re.findall("AF-(.*)-F1-model", x)[0]
        )

        # if it's the first results file, fill the dummy_df
        if i == 0:
            dummy_df = file_df
        # otherwise, add to the df
        else:
            dummy_df = pd.concat([dummy_df, file_df], axis=0)

    # extract unique uniprot IDs
    if dummy_df.empty:
        print(f"WARNING: No matching foldseek hits found in {input_files}.")
        hits = []
    else:
        hits = dummy_df["uniprotid"].unique()

    # if max_num_hits is provided, truncate the list
    if max_num_hits is not None:
        hits = hits[:max_num_hits]

    # save to a .txt file
    with open(output_file, "w+") as f:
        f.writelines(hit + "\n" for hit in hits)


# run this if called from the interpreter
def main():
    # parse arguments
    args = parse_args()

    # collect arguments individually
    input_files = args.input
    output_file = args.output
    evalue = args.evalue
    max_num_hits = args.max_num_hits

    # send to map_refseqids
    extract_foldseekhits(input_files, output_file, evalue=evalue, max_num_hits=max_num_hits)


# check if called from interpreter
if __name__ == "__main__":
    main()
