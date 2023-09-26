#!/usr/bin/env python
import argparse
import os

__all__ = ["prep_pdbpaths"]


# parse command line arguments
def parse_args():
    # Set command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--directories",
        nargs="+",
        help="Input path to directories to scrape PDB paths. Can accept multiple arguments.",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Name of output TXT file."
    )
    args = parser.parse_args()
    return args


def prep_pdbpaths(directories: list, output_file: str):
    """
    Takes a list of directories and converts it into a text file containing all PDB file paths.

    Args:
        directories (list): list of directory paths to scrape from
        output_file (str): path of destination file
    """
    output_text = []

    for directory in directories:
        pdbpaths = [
            os.path.join(directory, i) + "\n"
            for i in os.listdir(directory)
            if i.lower().endswith(".pdb")
        ]
        output_text = output_text + pdbpaths

    with open(output_file, "w+") as f:
        f.writelines(output_text)


def main():
    args = parse_args()
    directories = args.directories
    output_file = args.output

    prep_pdbpaths(directories, output_file)


if __name__ == "__main__":
    main()
