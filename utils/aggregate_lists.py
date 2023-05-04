#!/usr/bin/env python
import os
import argparse
from pathlib import Path

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs = '+', help = 'Paths of input files.')
    parser.add_argument("-o", "--output", help = 'Path of output file.')
    args = parser.parse_args()
    
    return args

# aggregate .txt files from folders based on suffix
def aggregate_lists(input_files: list, output_file: str):
    # empty set to collect ids and prevent collisions
    id_set = set()
    
    # iterate through files
    for file in input_files:
        # get file contents
        with open(file, 'r') as f:
            ids = f.read().splitlines()
            # add ids to set, which is non-redundant
            id_set.update(ids)
    
    # save unique entries to a new .txt file
    with open(output_file, 'w+') as f:
        f.writelines(id + '\n' for id in id_set)

# run this if called from the interpreter
def main():
    args = parse_args()
    input_files = args.input
    output_file = args.output
    aggregate_lists(input_files, output_file)

# check if called from interpreter
if __name__ == '__main__':
    main()
