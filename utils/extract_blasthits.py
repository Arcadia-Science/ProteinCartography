#!/usr/bin/env python3
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input")
parser.add_argument("-o", "--output")
args = parser.parse_args()

input_file = args.input
output_file = args.output

names = ['query_id', 'subject_id', 'per_identity', 'aln_length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'e-value', 'bit_score']
df = pd.read_csv(input_file, sep = '\t', names = names)

hits = df['subject_id'].unique()

with open(output_file, 'w+') as f:
    f.writelines(hit + '\n' for hit in hits)
    