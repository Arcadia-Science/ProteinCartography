#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np

# only import these functions when using import *
__all__ = ["extract_tmscore_feature"]

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = "Input distance matrix filepath.")
    parser.add_argument("-o", "--output", required = True, help = "Destination _features.tsv filepath.")
    parser.add_argument("-p", "--protid", required = True, help = "Protid of interest.")
    args = parser.parse_args()
    
    return args

def extract_tmscore_feature(input_file: str, protid: str, savefile = None) -> pd.DataFrame:
    '''
    Extracts the TMscore column from data for a specific protid and saves to file.
    This should be identical to the filename of the protein structure without the '.pdb' suffix.
    
    Args:
        input_file (str): path to input distance matrix.
        protid (str): protid value.
        savefile (str): path to destination file.
    '''
    # Read distance matrix
    df = pd.read_csv(input_file, sep = '\t', index_col = 'protid')
    
    # Get TMscore as features file (2-column: protid, TMscore_v_protid)
    df_dists = pd.DataFrame(df.loc[protid]).reset_index().rename(
        columns = {'index': 'protid', protid: f'TMscore_v_{protid}'})
    
    # Save to file if needed
    if savefile is not None:
        df_dists.to_csv(savefile, sep = '\t', index = None)
    
    # Return results
    return df_dists

# run this if called from the interpreter
def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    protid = args.protid
    
    extract_tmscore_feature(input_file, protid, savefile = output_file)
    
# check if called from interpreter
if __name__ == '__main__':
    main()
