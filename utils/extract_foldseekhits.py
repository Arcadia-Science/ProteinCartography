import pandas as pd
import re
import argparse
import os
from pathlib import Path

# default column names for a Foldseek run in this pipeline
FOLDSEEK_NAMES = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'prob', 'evalue', 'bits', 'qcov', 'tcov', 'qlan', 'taln', 'coord', 'tseq', 'taxid', 'taxname']

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs = '+', description = 'Takes .m8 file paths as input.')
    parser.add_argument("-o", "--output", description = 'Returns a .txt file as output.')
    args = parser.parse_args()
    
    return args

def extract_foldseekhits(input_files: str, output_file: str):
    # empty df for collecting results
    dummy_df = pd.DataFrame()
    
    # iterate through results files, reading them
    for i, file in enumerate(input_files):
        # load the file
        file_df = pd.read_csv(file, sep = '\t', names = FOLDSEEK_NAMES)
        
        # extract the model ID from the results target column
        file_df['modelid'] = file_df['target'].str.split(' ', expand = True)[0]
        
        # get the uniprot ID out from that target
        file_df['uniprotid'] = file_df['modelid'].apply(lambda x: re.findall('AF-(.*)-F1-model_v4', x)[0])
        
        # if it's the first results file, fill the dummy_df
        if i == 0:
            dummy_df = file_df
        # otherwise, add to the df
        else:
            dummy_df = pd.concat([dummy_df, file_df], axis = 0)
    
    # extract unique uniprot IDs
    hits = dummy_df['uniprotid'].unique()
    
    # save to a .txt file
    with open(output_file, 'w+') as f:
        f.writelines(hit + '\n' for hit in hits)

# run this if called from the interpreter
def main():
    # parse arguments
    args = parse_args()
    
    # collect arguments individually
    input_files = args.input
    output_file = args.output
    
    # send to map_refseqids
    extract_foldseekhits(input_files, output_file)

# check if called from interpreter
if __name__ == '__main__':
    main()
