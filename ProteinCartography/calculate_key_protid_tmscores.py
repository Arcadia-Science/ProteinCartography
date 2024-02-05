#!/usr/bin/env python
import argparse
import os
import subprocess
from pathlib import Path

from foldseek_clustering import pivot_foldseek_results

__all__ = [
    "run_foldseek_clustering",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-q",
        "--query-database",
        required=True,
        help="Path to the query database file.",
    )
    parser.add_argument(
        "-t",
        "--target-folder",
        required=True,
        help="Path to target folder that contains the target .pdb file.",
    )
    parser.add_argument(
        "-r",
        "--results-folder",
        required=True,
        help="Path to destination folder to save results.",
    )
    parser.add_argument(
        "-f",
        "--features-file",
        required=True,
        help="Path to the destination of the final TM-scores .csv file.",
    )
    args = parser.parse_args()

    return args


def run_foldseek_clustering(
    query_database: str,
    target_folder: str,
    results_folder: str,
    temp_folder=None,
    distances_filename="key_protid_tmscores.tsv",
):
    """
    Runs foldseek query_vs_target TMscore comparison of
    all query PDBs identified by the pipeline and a target PDB.
    Saves results to an output results_folder.
    Puts temporary files in a temp_folder
    (called `temp` in the results_folder if not specified explicity.)

    TODO (KC): Consider de-duplicating this method and `foldseek_clustering.run_foldseek_clustering`
    as it is nearly identical to this method.

    Args:
        query_database (str): path to the Foldseek database of the query .pdb files.
        target_folder (str): path to a folder with the target .pdb files.
        results_folder (str): path to a results folder.
        temp_folder (str): path to a temporary folder. Defaults to results_folder / temp_tm
        distances_filename (str): filename for output distances file.
            Defaults to `key_protid_tmscores.tsv`.
    Return:
        a tuple containing the full file paths of the distances file.
    """

    query_path = Path(query_database)
    target_path = Path(target_folder)
    results_path = Path(results_folder)

    foldseek_distances_tsv_query_vs_target = results_path / distances_filename

    if temp_folder is None:
        temp_path = results_path / "temp_tm"
    else:
        temp_path = Path(temp_folder)

    for path in [temp_path, results_path]:
        if not os.path.exists(path):
            os.mkdir(path)

    db_prefix_target = temp_path / "temp_db_target"
    subprocess.run(["foldseek", "createdb", target_path, db_prefix_target])

    foldseek_out_query_vs_target = temp_path / "query_vs_target"
    foldseek_tmp_query_vs_target = temp_path / "tmp_query_vs_target"
    subprocess.run(
        [
            "foldseek",
            "search",
            query_path,
            db_prefix_target,
            foldseek_out_query_vs_target,
            foldseek_tmp_query_vs_target,
            "-a",
            "--exhaustive-search",
        ]
    )

    foldseek_tmscore_query_vs_target = temp_path / "key_protid_tmscores.tsv"
    subprocess.run(
        [
            "foldseek",
            "aln2tmscore",
            query_path,
            db_prefix_target,
            foldseek_out_query_vs_target,
            foldseek_tmscore_query_vs_target,
        ]
    )

    subprocess.run(
        [
            "foldseek",
            "createtsv",
            query_path,
            db_prefix_target,
            foldseek_tmscore_query_vs_target,
            foldseek_distances_tsv_query_vs_target,
        ]
    )

    return str(foldseek_distances_tsv_query_vs_target)


def main():
    args = parse_args()
    query_database = args.query_database
    target_folder = args.target_folder
    results_folder = args.results_folder
    features_file = args.features_file

    # If there are no PDB files in the target_folder, then there is nothing for this script to do.
    # But to keep snakemake happy, we need to create the output file of the snakemake rule
    # that calls this script, and to keep downstream scripts happy,
    # the output file needs to be read by `pd.read_csv` as an empty dataframe.
    if not any(Path(target_folder).glob("*.pdb")):
        with open(features_file, "w") as file:
            file.write("protid\n")
            pass
        return

    distances_tsv = run_foldseek_clustering(query_database, target_folder, results_folder)
    pivot_foldseek_results(
        input_file=distances_tsv, output_file=features_file, column_prefix="TMscore_v_"
    )


if __name__ == "__main__":
    main()
