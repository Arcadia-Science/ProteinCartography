#!/usr/bin/env python
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = 'path of input list.')
    parser.add_argument("-o", "--output", required = True, help = 'path of output list.')
    args = parser.parse_args()
    
    return args

def rescue_mapping(input_file: str, output_file: str):
    with open(input_file, 'r') as f:
        text = [i.rstrip('\n') for i in f.readlines()]

    df = pd.DataFrame({'protid': text})
    unique = df['protid'].unique()

    with open(output_file, 'w+') as f:
        f.writelines([i + '\n' for i in unique])

def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    
    rescue_mapping(input_file, output_file)

if __name__ == '__main__':
    main()
