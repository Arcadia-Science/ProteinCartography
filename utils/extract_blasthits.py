#!/usr/bin/env python3
import argparse
import pandas as pd

FASTA_NAMES = ['query_id', 'subject_id', 'per_identity', 'aln_length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'e-value', 'bit_score']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    
    return args

# take an input blastresults file and create a .txt file from that
def extract_blasthits(input_file: str, output_file: str):
    df = pd.read_csv(input_file, sep = '\t', names = FASTA_NAMES)

    hits = df['subject_id'].unique()

    with open(output_file, 'w+') as f:
        f.writelines(hit + '\n' for hit in hits)
    
# run this if called from the interpreter
def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    extract_blasthits(input_file, output_file)

# check if called from interpreter
if __name__ == '__main__':
    main()
