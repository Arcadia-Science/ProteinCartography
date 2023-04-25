import pandas as pd
import re
import argparse
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input")
parser.add_argument("-o", "--output")
args = parser.parse_args()

input_directory = args.input
output_file = args.output

input_Path = Path(input_directory)
input_files = [input_Path / i for i in os.listdir(input_directory) if '.m8' in i]

dummy_df = pd.DataFrame()

names = ['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'prob', 'evalue', 'bits', 'qcov', 'tcov', 'qlan', 'taln', 'coord', 'tseq', 'taxid', 'taxname']
for i, file in enumerate(input_files):
    file_df = pd.read_csv(file, sep = '\t', names = names)
    file_df['modelid'] = file_df['target'].str.split(' ', expand = True)[0]
    file_df['uniprotid'] = file_df['modelid'].apply(lambda x: re.findall('AF-(.*)-F1-model_v4', x)[0])
    
    if i == 0:
        dummy_df = file_df
    else:
        dummy_df = pd.concat([dummy_df, file_df], axis = 0)

hits = dummy_df['uniprotid'].unique()
        
with open(output_file, 'w+') as f:
    f.writelines(hit + '\n' for hit in hits)
    