#!/usr/bin/env python
import pandas as pd
import re
import argparse
import os
from pathlib import Path

# only import these functions when using import *
__all__ = ["linear_convergence", "calculate_convergence"]

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tmscore-file", required = True, help = 'A TMscore-based distances TSV.')
    parser.add_argument("-f", "--fident-file", required = True, help = 'A fraction sequence identity-based distances TSV.')
    parser.add_argument("-p", "--protid", required = True, help = 'Unique identifier of the input protein.')
    parser.add_argument("-o", "--output", required = True, help = 'Returns a .tsv file of convergence features as output.')
    parser.add_argument("-m", "--method", default = 'linear', help = 'The method to use to calculate convergence/ divergence.')
    parser.add_argument("-e", "--evalue-maximum", default = '0.001', help = 'Maximum evalue cutoff. Default = 0.001')
    parser.add_argument("-T", "--tmscore-minimum", default = '0', help = 'Minimum tmscore cutoff. Default = 0')
    args = parser.parse_args()
    
    return args

def linear_convergence(tmscore, fident):
    distance = tmscore - fident
    return distance

def calculate_convergence(tmscore_file: str, fident_file: str, protid: str, output_file: str, 
                          method = 'linear', evalue_maximum = 0.001, tmscore_minimum = 0,
                          save = True):
    '''
    '''
    fident_df = pd.read_csv(fident_file, sep = '\t')
    tmscore_df = pd.read_csv(tmscore_file, sep = '\t')
    
    if method == 'linear':
        convergence_fxn = linear_convergence
    
    joint_df = fident_df.merge(tmscore_df, on = 'protid')
    joint_df[f'convergence_v_{protid}'] = joint_df.apply(lambda x: convergence_fxn(x[f'TMscore_v_{protid}'], x[f'fident_v_{protid}']), axis = 1)
    joint_df = joint_df[joint_df[f'evalue_v_{protid}'] < evalue_maximum]
    joint_df = joint_df[joint_df[f'TMscore_v_{protid}'] > tmscore_minimum]
    
    result_df = joint_df[['protid', f'convergence_v_{protid}']]
    
    if save:
        result_df.to_csv(output_file, index = None, sep = '\t')
    
    return result_df

# run this if called from the interpreter
def main():
    # parse arguments
    args = parse_args()
    
    # collect arguments individually
    tmscore_file = args.tmscore_file
    fident_file = args.fident_file
    protid = args.protid
    output_file = args.output
    method = args.method
    evalue_maximum = float(args.evalue_maximum)
    tmscore_minimum = float(args.tmscore_minimum)
    
    calculate_convergence(tmscore_file, fident_file, protid, output_file, 
                          method = method, evalue_maximum = evalue_maximum, tmscore_minimum = tmscore_minimum)

# check if called from interpreter
if __name__ == '__main__':
    main()
