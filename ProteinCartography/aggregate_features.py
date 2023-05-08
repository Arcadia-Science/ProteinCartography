#!/usr/bin/env python3
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs = '+', required = True, help = 'Paths of input features files.')
    parser.add_argument("-o", "--output", required = True, help = 'Path of output aggregated features file.')
    parser.add_argument("-v", "--override-file", default = '', help = 'Features override file for manual entries.')
    args = parser.parse_args()
    
    return args
    
def aggregate_features(input_files: list, output_file = '', features_override = ''):
    dfs = [pd.read_csv(file, sep = '\t') for file in input_files]
    
    agg_df = pd.DataFrame()
    
    for i, df in enumerate(dfs):
        if i == 0:
            agg_df = df
        else:
            agg_df = agg_df.merge(df, on = 'protid', how = 'outer')
    
    if features_override != '':
        features_override_df = pd.read_csv(features_override, sep = '\t')
        
        for entry in features_override_df['protid'].values:
            entry_row = features_override_df[features_override_df['protid'] == entry]
            if entry not in agg_df['protid'].values:
                agg_df = agg_df.add(entry_row)
            else:
                for col in [i for i in entry_row.columns if i != 'protid']:
                    if col not in agg_df.columns:
                        continue
                    
                    agg_df.loc[agg_df['protid'] == entry, col] = entry_row[col][0]
    
    if output_file != '':
        agg_df.to_csv(output_file, sep = '\t', index = None)
    
    return agg_df

def main():
    args = parse_args()
    input_files = args.input
    output_file = args.output
    features_override = args.override_file
    
    aggregate_features(input_files, output_file, features_override)
    
if __name__ == '__main__':
    main()
