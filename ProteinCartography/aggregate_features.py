#!/usr/bin/env python
import argparse

import pandas as pd

# only import these functions when using import *
__all__ = ["aggregate_features"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", nargs="+", required=True, help="Paths of input features files."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path of output aggregated features file."
    )
    parser.add_argument(
        "-v",
        "--override-file",
        default=None,
        help="Features override file for manual entries.",
        nargs="?",
    )
    args = parser.parse_args()

    return args


def aggregate_features(input_files: list, output_file=None, features_override=None) -> pd.DataFrame:
    """
    Aggregates features from multiple features files by `protid`.
    A features file contains at least one column with the name `protid`
    that is a unique identifier for each PDB file.
    The additional columns in the file are the features -
    the measures for each `protid` such as Length, LeidenCluster, etc.

    Returns a matrix where each of the columns is derived from the original features file
    and each row is a protein.

    Args:
        input_files (list): list of filepaths to aggregate.
        output_file (str): path of destination file.
        features_override (str): path to features override file.
            The override file contains as its first column a protid
            for specific proteins in the data.
            The remaining columns can be columns of any of the features files.
            After aggregation, the values of those specific cells will be updated
            to reflect the override file.
    Returns:
        a pandas.DataFrame containing the aggregated features.
    """

    # Read in list of dataframes
    dfs = [pd.read_csv(file, sep="\t") for file in input_files]

    # Create dummy dataframe for aggregation
    agg_df = pd.DataFrame()

    # Read the dataframes, merging them on protid
    for i, df in enumerate(dfs):
        if i == 0:
            agg_df = df
        else:
            agg_df = agg_df.merge(df, on="protid", how="outer")

    # If there's a features override file, use it
    if features_override is not None:
        # Read the file
        features_override_df = pd.read_csv(features_override, sep="\t")

        # Iterate over protids
        for entry in features_override_df["protid"].values:
            # Identify the row in the features override file where this protid is found
            entry_row = features_override_df[features_override_df["protid"] == entry]

            # If there isn't already a row for the data, append it
            if entry not in agg_df["protid"].values:
                # This should add np.nan for non-existing columns
                agg_df = agg_df.add(entry_row)
            # Otherwise, fill in the data at specific cells needed to be replaced
            else:
                # for each column in the override
                for col in [i for i in entry_row.columns if i != "protid"]:
                    # if the column doesn't exist, ignore it
                    if col not in agg_df.columns:
                        continue

                    # otherwise, replace the value at that position
                    agg_df.loc[agg_df["protid"] == entry, col] = entry_row[col].values[0]

    agg_df.drop_duplicates(inplace=True)

    # Save to an output file if a path is provided
    if output_file is not None:
        agg_df.to_csv(output_file, sep="\t", index=None)

    return agg_df


# run this if called from the interpreter
def main():
    args = parse_args()
    input_files = args.input
    output_file = args.output
    features_override = args.override_file

    aggregate_features(input_files, output_file, features_override)


# check if called from interpreter
if __name__ == "__main__":
    main()
