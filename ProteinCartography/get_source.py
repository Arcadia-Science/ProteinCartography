#!/usr/bin/env python
import argparse
import pandas as pd
import os

# only import these functions when using import *
__all__ = ["get_source"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input filepath containing list of protids as index. Usually pivoted distance matrix.",
    )
    parser.add_argument(
        "-f",
        "--hit-files",
        nargs="+",
        required=True,
        help="List of filepaths to parse source info.",
    )
    parser.add_argument("-o", "--output", required=True, help="Destination filepath.")
    parser.add_argument(
        "-k",
        "--keyids",
        nargs="+",
        default=[],
        help="Key ids (usually protids of input proteins) used to determine the source of a protid.\nUseful when there are multiple input proteins.",
    )
    args = parser.parse_args()

    return args


def get_source(
    input_file: str,
    hit_files: list,
    savefile=None,
    groupby=["method", "keyid"],
    methods=["blast", "foldseek"],
    keyids=[],
) -> pd.DataFrame:
    """
    Aggregates a list of hit files to determine the source of each protid for a features matrix.

    Args:
        input_file (str): path to input pivoted results file.
        hit_files (list): list of filepaths for results files (usually blasthits and foldseekhits).
        savefile (str): path to destination file.
        groupby (list): how to group summary metrics. Defaults to 'method' and 'keyid'.
        methods (list): what methods to look for Defaults to 'blast' and 'foldseek'.
        keyids (list): which keyids to identify as the source data for blast/foldseek.
    Returns:
        a pandas.DataFrame containing the results.
    """

    # Read the input file
    df = pd.read_csv(input_file, sep="\t", index_col="protid")

    # Get the index to get a list of protids
    df_indexes = pd.DataFrame(df.index)

    # Read the entries of each source file as a list
    for file in hit_files:
        # Get the name of the file
        sourcename = os.path.basename(file)
        # Get the column name (usually protid.method)
        sourcecol = sourcename.partition("hits")[0]

        # Read the file
        with open(file, "r") as f:
            sourceitems = [i.rstrip("\n") for i in f.readlines()]
        # For each protid, give it a 1 if it's a hit to that source file and a 0 if not.
        df_indexes[sourcecol] = pd.Series(
            [1 if i in sourceitems else 0 for i in df.index]
        )

    # Create a summary column by method.
    # For example, if there are multiple input proteins,
    # and you want to know if any of those produced a blast hit for each protid,
    # this will do that for you.
    if "method" in groupby:
        # Go across the methods
        for method in methods:
            # Identify the columns that have the method and see if any are true
            df_indexes[method] = df_indexes[df_indexes.filter(like=method).columns].any(
                axis=1
            )
            # convert the True / False to 1 / 0 encoding.
            df_indexes[method] = df_indexes[method].apply(
                lambda x: 1 if x == True else 0
            )
        # If both 'blast' and 'foldseek' are in the methods, generate a combined column for those hits
        # that were found via both methods.
        if "blast" in methods and "foldseek" in methods:
            df_indexes["blast+foldseek"] = df_indexes[["blast", "foldseek"]].all(axis=1)
            df_indexes["blast+foldseek"] = df_indexes["blast+foldseek"].apply(
                lambda x: 2 if x == True else 0
            )

        # Summarize across source columns to assign an "origin" for each protid.
        # If a protein came from both blast and foldseek,
        # because the blast+foldseek column has a max value of 2,
        # that assignment overrides the other two.
        df_indexes["source.method"] = df_indexes[
            ["blast", "foldseek", "blast+foldseek"]
        ].idxmax(axis=1)

    # Create a summary column by keyid (usually the input proteins used to make the search).
    # This tells you for each protid whether it was a hit to one of the keyid proteins.
    if "keyid" in groupby and keyids != []:
        for keyid in keyids:
            keyidcol = keyid + ".hit"
            df_indexes[keyidcol] = df_indexes[
                df_indexes.filter(like=keyid).columns
            ].any(axis=1)
            df_indexes[keyidcol] = df_indexes[keyidcol].apply(
                lambda x: 1 if x == True else 0
            )

    # Save to file
    if savefile is not None:
        df_indexes.to_csv(savefile, sep="\t", index=None)

    return df_indexes


# run this if called from the interpreter
def main():
    args = parse_args()
    input_file = args.input
    hit_files = args.hit_files
    output_file = args.output
    keyids = args.keyids

    get_source(input_file, hit_files, savefile=output_file, keyids=keyids)


# check if called from interpreter
if __name__ == "__main__":
    main()
