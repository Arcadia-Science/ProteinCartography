#!/usr/bin/env python
import argparse
import pandas as pd
import subprocess
import os
from pathlib import Path

# only import these functions when using import *
__all__ = [
    "run_foldseek_clustering",
    "make_struclusters_file",
    "pivot_foldseek_results",
]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-q",
        "--query-folder",
        required=True,
        help="Path to query folder containing all .pdb files of interest.",
    )
    parser.add_argument(
        "-r",
        "--results-folder",
        required=True,
        help="Path to destination folder to save results.",
    )
    args = parser.parse_args()

    return args


def run_foldseek_clustering(
    query_folder: str,
    results_folder: str,
    temp_folder=None,
    distances_filename="all_by_all_tmscore.tsv",
    cluster_filename="struclusters.tsv",
    cluster_mode="0",
    similarity_type="2",
):
    """
    Runs foldseek all-v-all TMscore comparison and clustering on all PDBs in an input query_folder.
    Saves results to an output results_folder.
    Puts temporary files in a temp_folder (called `temp` in the results_folder if not specified explicity.)

    Args:
        query_folder (str): path to a query folder containing .pdb files.
        results_folder (str): path to a results folder.
        temp_folder (str): path to a temporary folder. Defaults to results_folder / temp
        distances_filename (str): filename for output distances file. Defaults to `all_by_all_tmscore.tsv`.
        cluster_filename (str): filename for output struclusters file. Defaults to `struclusters.tsv`.
    Return:
        a tuple containing the full file paths of the distances file and the clusters file.
    """

    # Generate Path objects for each folder
    query_Path = Path(query_folder)
    results_Path = Path(results_folder)

    # count number of PDBs in query folder
    max_seqs = len([f for f in os.listdir(query_Path) if f.lower().endswith(".pdb")])

    # Generate Path object for temp folder
    if temp_folder is None:
        temp_Path = results_Path / "temp"
    else:
        temp_Path = Path(temp_folder)

    # Make Paths if they don't already exist
    for path in [temp_Path, results_Path]:
        if not os.path.exists(path):
            os.mkdir(path)

    # Run steps of Foldseek
    db_prefix = temp_Path / "temp_db"
    subprocess.run(["foldseek", "createdb", query_Path, db_prefix])

    foldseek_out = temp_Path / "all_by_all"
    foldseek_tmp = temp_Path / "tmp"
    subprocess.run(
        [
            "foldseek",
            "search",
            db_prefix,
            db_prefix,
            foldseek_out,
            foldseek_tmp,
            "-a",
            "max-seqs",
            str(max_seqs),
        ]
    )

    foldseek_tmscore = temp_Path / "all_by_all_tmscore"
    subprocess.run(
        [
            "foldseek",
            "aln2tmscore",
            db_prefix,
            db_prefix,
            foldseek_out,
            foldseek_tmscore,
        ]
    )

    foldseek_distancestsv = results_Path / distances_filename
    subprocess.run(
        [
            "foldseek",
            "createtsv",
            db_prefix,
            db_prefix,
            foldseek_tmscore,
            foldseek_distancestsv,
        ]
    )

    foldseek_cluster = temp_Path / "clu"
    subprocess.run(
        [
            "foldseek",
            "clust",
            db_prefix,
            foldseek_out,
            foldseek_cluster,
            "--cluster-mode",
            cluster_mode,
            "--similarity-type",
            similarity_type,
        ]
    )

    foldseek_clustertsv = results_Path / cluster_filename
    subprocess.run(
        [
            "foldseek",
            "createtsv",
            db_prefix,
            db_prefix,
            foldseek_cluster,
            foldseek_clustertsv,
        ]
    )

    # Return the output filepaths as a tuple
    return str(foldseek_distancestsv), str(foldseek_clustertsv)


def make_struclusters_file(foldseek_clustertsv: str, output_file: str):
    """
    Parses a Foldseek clusters file into a _features.tsv file.

    Args:
        foldseek_clustertsv (str): path of input clusters.tsv file.
        output_file (str): path of destination file.
    """
    # Read the input file
    df = pd.read_csv(foldseek_clustertsv, sep="\t", names=["ClusterRep", "protid"])

    # Strip the '.pdb' suffix so indices are protid
    df["ClusterRep"] = df["ClusterRep"].str.rstrip(".pdb")
    df["protid"] = df["protid"].str.rstrip(".pdb")

    # Aggregate groupings by ClusterRep
    df_merged = (
        df.groupby("ClusterRep")
        .agg(
            {
                i: ("first" if i == "ClusterRep" else lambda x: [i for i in x])
                for i in df.columns
            }
        )
        .reset_index(drop=True)
    )

    # Determine max number of characters for padding purposes
    max_chars = len(str(df_merged.index[-1]))

    # Generate structural cluster IDs and make a new column with those IDs
    SC_ids = "SC" + pd.Series(df_merged.index).apply(lambda x: str(x).zfill(max_chars))
    df_merged.insert(0, "StruCluster", SC_ids)

    # Drop the representative cluster member column
    df_merged.drop(columns=["ClusterRep"], inplace=True)

    # Explode the data into a one-protid-per-rown dataframe
    df_exploded = df_merged.explode("protid")

    # Slice only the relevant columns and save to a tsv file.
    df_exploded = df_exploded[["protid", "StruCluster"]]
    df_exploded.to_csv(output_file, sep="\t", index=None)

    return df_exploded


def pivot_foldseek_results(input_file: str, output_file: str):
    """
    Takes a df containing cleaned foldseek results and creates a similarity matrix.

    Args:
        input_file (str): input cleaned foldseek results filepath
        output_file (str): output similarity matrix filepath
    """
    import pandas as pd
    import os

    # read the foldseek output with delimited whitespace
    tmscore_df = pd.read_csv(input_file, delim_whitespace=True, header=None)
    tmscore_df = tmscore_df[[0, 1, 2]]
    foldseek_df = tmscore_df.rename(columns={0: "protid", 1: "target", 2: "tmscore"})
    foldseek_df.drop_duplicates(["protid", "target"], inplace=True)

    # pivot the data so that it's a square matrix, filling empty comparisons with 0
    pivoted_table = pd.pivot(
        foldseek_df, index="protid", columns="target", values="tmscore"
    ).fillna(0)

    # remove .pdb from row and column indices
    pivoted_table.columns = pivoted_table.columns.str.removesuffix(".pdb")
    pivoted_table.index = pivoted_table.index.str.removesuffix(".pdb")

    # save to file
    pivoted_table.to_csv(output_file, sep="\t")

    return pivoted_table


# run this if called from the interpreter
def main():
    args = parse_args()
    query_folder = args.query_folder
    results_folder = args.results_folder

    distancestsv, clusterstsv = run_foldseek_clustering(query_folder, results_folder)

    pivotedtsv = distancestsv.replace(".tsv", "_pivoted.tsv")
    pivot_foldseek_results(distancestsv, pivotedtsv)

    featurestsv = clusterstsv.replace(".tsv", "_features.tsv")
    make_struclusters_file(clusterstsv, featurestsv)


# check if called from interpreter
if __name__ == "__main__":
    main()
