#!/usr/bin/env python
import argparse

# only import these functions when using import *
__all__ = ["aggregate_lists"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        nargs="+",
        help="Paths of input files. Each file should be a .txt file with one accession per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path of output file. Duplicates across the lists will be removed.",
    )
    args = parser.parse_args()
    return args


# aggregate .txt files from folders based on suffix
def aggregate_lists(input_files: list, output_file: str):
    """
    Takes in a list of input .txt files containing accessions and combines them,
    removing duplicates.

    Args:
        input_files (list): list of string paths of input files.
        output_file (str): path of destination file.
    """
    # empty set to collect ids and prevent collisions
    id_set = set()

    # iterate through files
    for file in input_files:
        # get file contents
        with open(file) as f:
            ids = f.read().splitlines()
            # add ids to set, which is non-redundant
            id_set.update(ids)

    # save unique entries to a new .txt file
    with open(output_file, "w+") as f:
        f.writelines(id + "\n" for id in id_set)


# run this if called from the interpreter
def main():
    args = parse_args()
    input_files = args.input
    output_file = args.output
    aggregate_lists(input_files, output_file)


# check if called from interpreter
if __name__ == "__main__":
    main()
