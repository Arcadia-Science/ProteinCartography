#!/usr/bin/env python
import os
import argparse
from pathlib import Path

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs = '+', help = 'Names of input directories.')
    parser.add_argument("-s", "--suffix", nargs = '+', help = 'Suffix of each output file in each directory.')
    parser.add_argument("-o", "--output", help = 'Name of output file.')
    args = parser.parse_args()
    
    return args

# aggregate .txt files from folders based on suffix
def aggregate_lists(input_dirs: list, input_suffixes: list, output_file: str):
    
    # make sure that input directories and suffixes are of same length
    if len(input_dirs) != len(input_suffixes):
        sys.exit('Number of input directories and suffixes does not match. Please try again.')
        
    # zip input directories and suffixes into dict
    input_dict = dict(zip(input_dirs, input_suffixes))

    # empty list to collect file paths
    files = []
    
    # iterate through input dictionary
    for directory, suffix in input_dict.items():
        # make Path object for location
        directory_Path = Path(directory)
        # generate list of filepaths
        matches = [directory_Path / i for i in os.listdir(directory) if suffix in i]
        # append to files list
        files = files + matches

    # empty set to collect ids and prevent collisions
    id_set = set()
    
    # iterate through files
    for file in files:
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
    input_dirs = args.input
    input_suffixes = args.suffix
    output_file = args.output
    aggregate_lists(input_dirs, input_suffixes, output_file)

# check if called from interpreter
if __name__ == '__main__':
    main()
    