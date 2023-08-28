#!/usr/bin/env python
import pandas as pd
import re
import argparse
import os
from pathlib import Path

# only import these functions when using import *
__all__ = ["aggregate_foldseek_fident"]

# default column names for a Foldseek run in this pipeline
FOLDSEEK_NAMES = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'prob', 'evalue', 'bits', 'qcov', 'tcov', 'qlan', 'taln', 'coord', 'tseq', 'taxid', 'taxname']

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs = '+', required = True, help = 'Takes .m8 file paths as input.')
    parser.add_argument("-o", "--output", required = True, help = 'Returns a .tsv file as output.')
    parser.add_argument("-p", "--protid", required = True, help = 'Unique protein identifier used to generate .m8 files from search.')
    args = parser.parse_args()
    
    return args

def aggregate_foldseek_fident(input_files: list, output_file: str, protid: str) -> pd.DataFrame:
    '''
    Takes a list of input tabular Foldseek results files from the API query (ending in .m8).
    Generates a fident_features.tsv file containing the fraction sequence identity returned by Foldseek.
    
    Args:
        input_files (list): list of string paths to input files.
        output_file (str): path of destination file.
        protid (str): protid associated with the input files.
    '''
    
    # Create a collector dataframe
    dummy_df = pd.DataFrame()
    
    # open each .m8 file and append it to the collector, if it's not empty
    for f in input_files:
        temp_df = pd.read_csv(f, sep = '\t', names = FOLDSEEK_NAMES)
        
        if len(temp_df) == 0:
            continue
        
        dummy_df = pd.concat([dummy_df, temp_df])
    
    # tidy up index and any duplicates
    dummy_df.reset_index(inplace = True, drop = True)
    dummy_df.drop_duplicates(inplace = True)

    # extract the model ID from the results file
    dummy_df['modelid'] = dummy_df['target'].str.split(' ', expand = True)[0]
    dummy_df = dummy_df[dummy_df['modelid'].str.contains('-F1-model_v4')]
    
    # get the uniprot ID out from that target
    dummy_df['protid'] = dummy_df['modelid'].apply(lambda x: re.findall('AF-(.*)-F1-model_v4', x)[0])
    
    # collect only the relevant columns
    results_df = dummy_df[['protid', 'fident', 'prob', 'evalue']]
    results_df = results_df.drop_duplicates()
    
    # rescale fraction identity from percentile to actual fraction
    results_df[f'fident_v_{protid}'] = results_df['fident'] * 0.01
    results_df.drop(columns = 'fident', inplace = True)
    
    # append protid to evalue and prob columns
    results_df.rename(columns = {'evalue': f'evalue_v_{protid}', 'prob': f'prob_v_{protid}'}, inplace = True)
    
    # save to TSV
    results_df.to_csv(output_file, sep = '\t', index = None)
    
    # return Pandas DataFrame
    return results_df

# run this if called from the interpreter
def main():
    # parse arguments
    args = parse_args()
    
    # collect arguments individually
    input_files = args.input
    output_file = args.output
    protid = args.protid
    
    # send to aggregate_foldseek_fident
    aggregate_foldseek_fident(input_files, output_file, protid)

# check if called from interpreter
if __name__ == '__main__':
    main()
