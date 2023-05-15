#!/usr/bin/env python
import argparse
import pandas as pd
from bioservices import UniProt

# only import these functions when using import *
__all__ = ["query_uniprot_bioservices"]

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = "path of input merged hits file")
    parser.add_argument("-o", "--output", required = True, help = "path of destination uniprot_features.tsv file")
    args = parser.parse_args()
    
    return args

def query_uniprot_bioservices(input_file: str, output_file: str, save = True) -> pd.DataFrame:
    '''
    Takes an input list of accessions and gets the full information set from Uniprot for those proteins.
    
    Args:
        input_file (str): path of input list text file where each accession is on a new line
        output_file (str): path of destination tsv file with all uniprot features
    
    Returns:
        a pandas.DataFrame of the resulting features
    '''
    
    # open and read the file
    with open(input_file, 'r') as f:
        id_list = [i.rstrip('\n') for i in f.readlines()]
    
    # perform ID mapping using bioservices UniProt
    # should probably do this differently in the future because it often results in weird memory leaks
    u = UniProt()
    results = u.mapping("UniProtKB_AC-ID", "UniProtKB", query = ' '.join(id_list))
    
    # read the results as a normalized json
    results_df = pd.json_normalize(results['results'])
    # remove the "to" prefix for tidier columns
    results_df.columns = results_df.columns.str.removeprefix('to.')
    # add a protid column for later merging
    results_df.insert(0, 'protid', results_df['primaryAccession'])
    
    # save to file if needed
    if save:
        results_df.to_csv(output_file, index = None, sep = '\t')
    
    return results_df

# run this if called from the interpreter
def main():
    args = parse_args()
    
    input_file = args.input
    output_file = args.output
    
    query_uniprot_bioservices(input_file, output_file)

# check if called from interpreter
if __name__ == '__main__':
    main()
