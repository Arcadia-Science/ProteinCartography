import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True)
    parser.add_argument("-o", "--output", required = True)
    parser.add_argument("-p", "--protid", required = True)
    args = parser.parse_args()
    
    return args

def extract_tmscore_feature(input_file: str, protid: str, savefile = ''):
    df = pd.read_csv(input_file, sep = '\t', index_col = 'protid')
    df_dists = pd.DataFrame(df.loc[protid]).reset_index().rename(
        columns = {'index': 'protid', protid: f'TMscore_v_{protid}'})
    
    if savefile != '':
        df_dists.to_csv(savefile, sep = '\t', index = None)
    
    return df_dists

def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    protid = args.protid
    
    extract_tmscore_feature(input_file, protid, savefile = output_file)
    
if __name__ == '__main__':
    main()
