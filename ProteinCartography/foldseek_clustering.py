#!/usr/bin/env python
import argparse
from collections import defaultdict
import csv
import pandas as pd
import subprocess
import os
from pathlib import Path

# only import these functions when using import *
__all__ = [
    "run_foldseek_clustering",
    "make_struclusters_file",
    "reading_data",
    "get_line_for_protid",
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
        ["foldseek", "search", db_prefix, db_prefix, foldseek_out, foldseek_tmp, "-a"]
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


def reading_data(input_file: str):
    """
    Read in a cleaned foldseek results file.
    Return a dictionary whose keys are protids and whose entries
    are another dictionary of all the targets to their corresponding
    tmscore. This partitioning by protid allows for easier parallelization.

    Args:
        input_file (str): input cleaned foldseek results filepath
    Return:
        A tuple containing protid to target mappings and all targets in dataset.
    """
    targets = set()
    entries = defaultdict(dict)

    with open(input_file) as fh:
        for line in fh:
            # Getting the first three entries in the output from Foldseek using
            # the parameters specified in `run_foldseek_clustering`, the full
            # list of params avail is specified https://github.com/steineggerlab/foldseek#output-search
            protid, target, score, *_ = [e.strip() for e in line.split()]
            protid = protid.replace(".pdb", "")
            target = target.replace(".pdb", "")

            if target in entries[protid]:
                if entries[protid][target] != score:
                    raise ValueError(
                        f"Multiple values supplied for protid={protid}, target={target} with different scores."
                    )
            else:
                entries[protid][target] = score
            targets.add(target)

    return entries, targets


def get_line_for_protid(protid_and_targets: tuple, targets: set):
    """
    Given a protid_and_targets tuple and a list of all possible targets, this function returns
    the line of the similarity matrix file for the given protid in the entry.

    Args:
        protid_and_targets (tuple): the first entry is a protid and the second entry is a target to score dictionary.
        targets (set): A set of all targets seen in the data set. This allows setting 0.0 as fillna(0.0) did with pandas.
    Return:
        A list with the protid as the first element and then scores for each target.
    """
    protid, targets_to_scores = protid_and_targets
    scores = []
    for target in targets:
        scores.append(targets_to_scores.get(target, "0.0"))
    return [protid] + scores


def pivot_foldseek_results(input_file: str, output_file: str):
    """
    Takes a file with the first three columns being protid, target, and the tmscore.
    It then saves a similarity matrix to a csv. There is no return value.

    Args:
        input_file (str): input cleaned foldseek results filepath
        output_file (str): output similarity matrix filepath
    Return:
        None
    """
    entries, targets = reading_data(input_file)

    with open(output_file, "w", newline="") as fh:
        csv_writer = csv.writer(fh, delimiter="\t")

        header = ["protid"] + list(targets)
        csv_writer.writerow(header)

        for entry in sorted(entries.items()):
            csv_writer.writerow(get_line_for_protid(entry, targets))


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
