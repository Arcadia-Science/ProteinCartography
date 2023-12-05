#!/usr/bin/env python
import argparse

# only import these functions when using import *
__all__ = ["rescue_mapping"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="path of input list.")
    parser.add_argument("-o", "--output", required=True, help="path of output list.")
    args = parser.parse_args()

    return args


def rescue_mapping(input_file: str, output_file: str):
    """
    Takes an input file with accessions that may not be unique, and makes them unique.

    Args:
        input_file (str): path of input text file with one accession per line
        output_file (str): path of destination file to write sanitized accessions
    """

    # read the file
    with open(input_file) as f:
        text = [i.rstrip("\n") for i in f.readlines()]

    # unique-ify it
    unique = set(text)

    # save that set
    with open(output_file, "w+") as f:
        f.writelines([i + "\n" for i in unique])


# run this if called from the interpreter
def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output

    rescue_mapping(input_file, output_file)


# check if called from interpreter
if __name__ == "__main__":
    main()
