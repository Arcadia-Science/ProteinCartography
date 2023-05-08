import argparse
import subprocess
import os
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, nargs = '+', help = 'Input file path of a .txt file with one accession per line.')
    parser.add_argument("-o", "--output", required = True, help = 'Output directory to save dummy files into.')
    parser.add_argument("-M", "--maximum", default = '0', help = 'Maximum number of dummies to generate.\nThis will cause Snakemake to download only this many entries.')
    args = parser.parse_args()
    return args

def make_dummies(input_file: str, output_dir: str, maximum = 0):
    '''
    Creates empty files matching the Alphafold accessions for Snakemake to create wildcards with.
    
    Args:
        input_file (str): path to input .txt file with one accession per line.
        output_dir (str): path to output directory.
        maximum (int): maximum number of entries to create dummy files for. If 0, creates all.
    '''
    if output_dir[-1] != '/':
        output_dir = output_dir + '/'
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    input_text = []
    with open(input_file[0], 'r') as f:
        input_text = f.read().splitlines()
    
    if maximum == 0 or maximum > len(input_text):
        pass
    else:
        input_text = input_text[:maximum]
    
    for acc in input_text:
        output_path = output_dir + acc + '.txt'
        Path(output_path).touch()

def main():
    args = parse_args()
    input_file = args.input
    output_dir = args.output
    maximum = int(args.maximum)
    make_dummies(input_file, output_dir, maximum)
    
if __name__ == '__main__':
    main()
