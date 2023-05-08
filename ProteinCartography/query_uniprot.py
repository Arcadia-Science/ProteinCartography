#!/usr/bin/env python3
import argparse
import pandas as pd
from bioservices import UniProt

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True)
    parser.add_argument("-o", "--output", required = True)
    args = parser.parse_args()
    
    return args

def query_uniprot(input_file: str, output_file: str, save = True):
    with open(input_file, 'r') as f:
        id_list = [i.rstrip('\n') for i in f.readlines()]
    
    u = UniProt()
    results = u.mapping("UniProtKB_AC-ID", "UniProtKB", query = ' '.join(id_list))
    
    results_df = pd.json_normalize(results['results'])
    results_df.columns = results_df.columns.str.removeprefix('to.')
    results_df.insert(0, 'protid', results_df['primaryAccession'])
    
    if save:
        results_df.to_csv(output_file, index = None, sep = '\t')
    
    return results_df

def main():
    args = parse_args()
    
    input_file = args.input
    output_file = args.output
    
    query_uniprot(input_file, output_file)

if __name__ == '__main__':
    main()
