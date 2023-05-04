#!/usr/bin/env python
from bioservices import UniProt
import argparse
import pandas as pd

# check through these default databases
DEFAULT_DBS = ['EMBL-GenBank-DDBJ_CDS', 'RefSeq_Protein']

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-d", "--databases", nargs = '+', default = DEFAULT_DBS)
    args = parser.parse_args()
    return args

# takes a list of IDs and maps them to Uniprot using bioservices
# might make a more generalizable version of this and put it somewhere else
def map_refseqids(input_file: str, output_file: str, query_dbs: list):
    
    # make object that references UniProt database
    u = UniProt()
    
    # open the input file to extract ids
    with open(input_file, 'r') as f:
        ids = f.read().splitlines()
        
    if len(ids) > 100000:
        ids = ids[0:100000]
        
    # make an empty collector dataframe for mapping
    dummy_df = pd.DataFrame()
    
    # for each query database, map
    for i, db in enumerate(query_dbs):
        # u.mapping returns a gross json file
        results = u.mapping(db, "UniProtKB", query = ' '.join(ids))
        
        # pandas can normalize the json and make it more tractable
        results_df = pd.json_normalize(results['results'])
        
        # if there are no results, move on
        if len(results_df) == 0:
            continue
        
        # if it's the first database, replace it with the dummy dataframe
        if i == 0:
            dummy_df = results_df
        # otherwise append to the dataframe
        else:
            dummy_df = pd.concat([dummy_df, results_df], axis = 0)
    
    # extract just the unique Uniprot accessions
    hits = dummy_df['to.primaryAccession'].unique()
    
    # save those accessions to a .txt file
    with open(output_file, 'w+') as f:
        f.writelines(hit + '\n' for hit in hits)

# run this if called from the interpreter
def main():
    # parse arguments
    args = parse_args()
    
    # collect arguments individually
    input_file = args.input
    output_file = args.output
    query_dbs = args.databases
    
    # send to map_refseqids
    map_refseqids(input_file, output_file, query_dbs)

# check if called from interpreter
if __name__ == '__main__':
    main()
