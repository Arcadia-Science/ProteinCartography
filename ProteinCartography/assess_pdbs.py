#!/usr/bin/env python
import argparse
import re
import pandas as pd
import numpy as np
import arcadia_pycolor as apc
from io import StringIO
import os

__all__ = [
    'fetch_atoms', 'fetch_dbref', 
    'fetch_experiment', 'fetch_title', 
    'extract_residue_confidence', 'assign_residue_colors',
    'parse_chains', 'assign_origin', 'read_txtlist',
    'assess_pdbs'
]

RESIDUE_CONFIDENCE_COLORS = {
    'very_high': '#4A72B0',
    'confident': apc.All['arcadia:vitalblue'],
    'low': apc.All['arcadia:canary'],
    'very_low': apc.All['arcadia:amber']
}

RESIDUE_BINS = ['very_low', 'low', 'confident', 'very_high']

# based on: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
ATOM_SPEC_DICT = {
    'ATOM': (0, 6), 
    'SERIAL': (6, 11), 
    'NAME': (12, 16), 
    'ALTLOC': (16, 17), 
    'RESIDUE': (17, 20), 
    'CHAIN': (21, 22), 
    'RESNUM': (22, 26), 
    'RESINS': (26, 27), 
    'X': (30, 38), 
    'Y': (38, 45), 
    'Z': (46, 54), 
    'OCC': (54, 60), 
    'TEMP': (60, 66), 
    'SEG': (72, 76), 
    'ELEM': (76, 78), 
    'CHARGE': (78, 80)
}

ATOM_SPEC_VALS = list(ATOM_SPEC_DICT.values())
ATOM_SPEC_NAMES = list(ATOM_SPEC_DICT.keys())

# parse command line arguments
def parse_args():
    # Set command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--textfile", help = 'Input path to text file containing file paths, one per line.')
    parser.add_argument("-o", "--output", required = True, help = 'Name of output TSV file.')
    args = parser.parse_args()
    return args

def fetch_atoms(input_path: str) -> pd.DataFrame:
    
    with open(input_path, 'r') as f:
        atoms = [i for i in f.readlines() if 'ATOM' in i[0:6]]
    
    if len(atoms) == 0:
        return pd.DataFrame()
        
    data = pd.read_fwf(StringIO(''.join(atoms)), names = ATOM_SPEC_NAMES, colspecs = ATOM_SPEC_VALS)
    
    return data

def fetch_dbref(input_path: str) -> pd.DataFrame:
    with open(input_path, 'r') as f:
        dbref = [i for i in f.readlines() if 'DBREF' in i[0:6]]
        
    if len(dbref) == 0:
        return pd.DataFrame()
        
    data = pd.read_fwf(StringIO(''.join(dbref)), header = None)
    
    return data

def fetch_experiment(input_path: str):
    with open(input_path, 'r') as f:
        expdta = [' '.join([i for i in i.split() if i != 'EXPDTA']) for i in f.readlines() if 'EXPDTA' in i]
        
    expdta_out = ' '.join(expdta).split(';')
    
    return expdta_out

def fetch_title(input_path: str):
    with open(input_path, 'r') as f:
        title = [' '.join([i for i in i.split() if i != 'TITLE']) for i in f.readlines() if 'TITLE' in i]
        
    title_out = ' '.join(title)
    
    return title_out

def extract_residue_confidence(input_path: str):
    data = fetch_atoms(input_path)
    
    return list(data['TEMP'].astype(float))

def assign_residue_colors(lst: list):
    bins = [0, 50, 70, 90, 100]
    result = np.digitize(lst, bins, right=True)
    colors = [RESIDUE_CONFIDENCE_COLORS[RESIDUE_BINS[i - 1]] for i in result]
    
    return colors

def parse_chains(input_path: str):
    data = fetch_atoms(input_path)
    chains = list(data['CHAIN'].unique())
    
    return chains

def assign_origin(input_path: str):
    
    AF_FLAG, AF_TITLE_FLAG = 0, 0
    PDB_FLAG, PDB_REF_FLAG = 0, 0
    ESM_FLAG, ESM_TITLE_FLAG = 0, 0
    
    with open(input_path, 'r') as f:
        contents = f.read()
        
        if 'ALPHAFOLD' in contents:
            AF_FLAG = 1
        if 'PDB' in contents:
            PDB_FLAG = 1
        if 'ESMFOLD' in contents:
            ESM_FLAG = 1
    try:
        dbref = fetch_dbref(input_path)[5].values
    except KeyError:
        print(f'{input_path} failed')
        dbref = []
    
    if 'PDB' in dbref:
        PDB_REF_FLAG = 2
    
    title = fetch_title(input_path)
    if 'ALPHAFOLD' in title:
        AF_TITLE_FLAG = 2
    elif 'ESMFOLD' in title:
        ESM_TITLE_FLAG = 2
    
    AF_SCORE = AF_FLAG + AF_TITLE_FLAG
    PDB_SCORE = PDB_FLAG + PDB_REF_FLAG
    ESM_SCORE = ESM_FLAG + ESM_TITLE_FLAG
    
    OTHER_SCORE = 3 - AF_SCORE - PDB_SCORE - ESM_SCORE
    
    scores = {'Alphafold': AF_SCORE, 'ESMFold': ESM_SCORE, 'PDB': PDB_SCORE, 'Other': OTHER_SCORE}
    
    maxscore = max(zip(scores.values(), scores.keys()))[1]
    
    return maxscore

def assess_pdbs(structures_list: list, output_file = None):
    collector_df = pd.DataFrame()
    
    for structure_file in structures_list:
        if not os.path.exists(structure_file):
            continue
        
        with open(structure_file, 'r') as f:
            if '<Error>' in f.read():
                continue
        
        expdta = fetch_experiment(structure_file)
        
        origin = assign_origin(structure_file)
        
        if origin != 'PDB':
            confidence = np.mean(extract_residue_confidence(structure_file))
        else:
            confidence = 100
            
        chains = parse_chains(structure_file)
            
        info_df = pd.DataFrame({
            'protid': os.path.basename(structure_file).split('.pdb')[0],
            'pdb_origin': [origin], 
            'pdb_confidence': [confidence], 
            'pdb_chains': [chains],
        })
        
        collector_df = pd.concat([collector_df, info_df], ignore_index = True)
        
    if output_file is not None:
        collector_df.to_csv(output_file, sep = '\t', index = None)
        
    return collector_df

def read_txtlist(input_file: str):
    with open(input_file, 'r') as f:
        paths = [i.rstrip('\n') for i in f.readlines()]
        
    return paths

def main():
    args = parse_args()
    textfile = args.textfile
    output_file = args.output
    
    files_list = read_txtlist(textfile)
    
    assess_pdbs(structures_list = files_list, output_file = output_file)
    
if __name__ == '__main__':
    main()
