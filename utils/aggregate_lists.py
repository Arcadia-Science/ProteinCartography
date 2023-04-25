#!/usr/bin/env python3
import os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", nargs = '+', help = 'Names of input directories.')
parser.add_argument("-s", "--suffix", nargs = '+', help = 'Suffix of each output file in each directory.')
parser.add_argument("-o", "--output", help = 'Name of output file.')
args = parser.parse_args()

output_file = args.output

if len(input_dirs := args.input) != len(input_suffixes := args.suffix):
    sys.exit('Number of input directories and suffixes does not match. Please try again.')

input_dict = dict(zip(input_dirs, input_suffixes))

files = []
for directory, suffix in input_dict.items():
    directory_Path = Path(directory)
    matches = [directory_Path / i for i in os.listdir(directory) if suffix in i]
    files = files + matches

id_set = set()

for file in files:
    with open(file, 'r') as f:
        ids = f.read().splitlines()
        id_set.update(ids)

with open(output_file, 'w+') as f:
    f.writelines(id + '\n' for id in id_set)
    