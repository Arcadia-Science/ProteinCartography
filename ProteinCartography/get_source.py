import argparse
import pandas as pd
import numpy as np
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True)
    parser.add_argument("-f", "--hit-files", nargs = '+', required = True)
    parser.add_argument("-o", "--output", required = True)
    parser.add_argument("-k", "--keyids", nargs = '+', default = [])
    args = parser.parse_args()
    
    return args

def get_source(input_file: str, hit_files: list, 
               savefile = '', groupby = ['method', 'keyid'], 
               methods = ['blast', 'foldseek'], keyids = []):
    df = pd.read_csv(input_file, sep = '\t', index_col = 'protid')
    df_indexes = pd.DataFrame(df.index)
    
    for file in hit_files:
        sourcename = os.path.basename(file)
        sourcecol = sourcename.partition('hits')[0]
        with open(file, 'r') as f:
            sourceitems = [i.rstrip('\n') for i in f.readlines()]
        df_indexes[sourcecol] = pd.Series([1 if i in sourceitems else 0 for i in df.index])
        
    if 'method' in groupby:
        for method in methods:
            df_indexes[method] = df_indexes[df_indexes.filter(like=method).columns].any(axis = 1)
            df_indexes[method] = df_indexes[method].apply(lambda x: 1 if x == True else 0)
        if 'blast' in methods and 'foldseek' in methods:
            df_indexes['blast+foldseek'] = df_indexes[['blast', 'foldseek']].all(axis = 1)
            df_indexes['blast+foldseek'] = df_indexes['blast+foldseek'].apply(lambda x: 2 if x == True else 0)
        
        df_indexes['source.method'] = df_indexes[['blast', 'foldseek', 'blast+foldseek']].idxmax(axis = 1)
    
    if 'keyid' in groupby and keyids != []:
        for keyid in keyids:
            keyidcol = keyid + '.hit'
            df_indexes[keyidcol] = df_indexes[df_indexes.filter(like=keyid).columns].any(axis = 1)
            df_indexes[keyidcol] = df_indexes[keyidcol].apply(lambda x: 1 if x == True else 0)
            
    if savefile != '':
        df_indexes.to_csv(savefile, sep = '\t', index = None)

    return df_indexes

def main():
    args = parse_args()
    input_file = args.input
    hit_files = args.hit_files
    output_file = args.output
    keyids = args.keyids
    
    get_source(input_file, hit_files, savefile = output_file, keyids = keyids)
    
if __name__ == '__main__':
    main()
    