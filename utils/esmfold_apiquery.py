#!/usr/bin/env python3
import sys
import os
import argparse
from requests import post
from time import sleep

### NOTES
#ESMFold API example from website:
"""
curl -X POST --data "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL" https://api.esmatlas.com/foldSequence/v1/pdb/
"""

FASTA_FORMATS = ['fa', 'fna', 'fasta']

def parse_args():
    # Set command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = 'Name of input file. Must be a single-entry peptide FASTA file (ends with .fa, .fna, .fasta).')
    parser.add_argument("-o", "--output", default = '',
                        help = 'Name of output file. If not provided, replaces ".fasta" from input with ".pdb".')
    args = parser.parse_args()
    return args

def post_query(fasta: str):
    result = post('https://api.esmatlas.com/foldSequence/v1/pdb/',
                  data = fasta)
    if result.status_code == 200:
        return result.text
    else:
        print(f"Error: {result.status_code}")
        return None

def esmfold_apiquery(input_file: str, output_file = ''):
    # Check to make sure input file has a correct FASTA suffix
    if not any([fmt in input_file.rsplit('.', 1)[1].lower() for fmt in FASTA_FORMATS]):
        sys.exit(f'Input expects a FASTA file ({FASTA_FORMATS})')

    # Makes sure that the input file exists
    if not os.path.exists(input_file):
        sys.exit(f'File {input_file} not found.')
        
    # if output file is not provided, generate a name for output file
    if output_file == '':
        output_filepath = input_file.replace(input_file.rsplit('.', 1)[1], 'pdb')
    else:
        output_filepath = output_file

    # Open input file and collect text as string
    with open(input_file, 'r') as f:
        content = f.read()
        
        entries = content.split(">")[1:]
        if len(entries) > 1:
            raise Exception('This script expects a single FASTA entry in the input file. Please try again.')
        
        lines = entries[0].strip().split("\n")
        fasta = "".join(lines[1:])
    
    if (prot_len := len(fasta)) > 400:
        print(f'The input protein is {prot_len} AA long.\nESMFold API query only allows proteins up to 400 AA long.\nTry using ColabFold instead.\nSkipping...')
        return
    
    # submit a new job via the API
    result = post_query(fasta)
    
    if result is not None:
        with open(output_filepath, 'w+') as file:
            file.write(result)

# run this if called from the interpreter
def main():
    # parse args
    args = parse_args()
    
    input_file = args.input
    output_file = args.output
    
    esmfold_apiquery(input_file, output_file)

# check if called from interpreter
if __name__ == '__main__':
    main()
