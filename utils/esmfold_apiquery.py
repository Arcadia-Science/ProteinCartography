#!/usr/bin/env python3
import sys
import os
import argparse
from requests import get, post
from time import sleep

FASTA_FORMATS = ['fa', 'fna', 'fasta']

def parse_args():
    # Set command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = 'Name of input file. Must be a FASTA file.')
    parser.add_argument("-o", "--output", default = '',
                        help = 'Name of output file. If not provided, replaces ".fasta" from input with ".pdb".')
    args = parser.parse_args()
    return args

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

    # Collector for PDB information for requests.post()
    fasta = ''

    # Open input file and collect text as string
    with open(input_file, 'r') as file:
        text = file.readlines()
        
        if len([line for line in text if '>' in line]) > 1:
            raise Exception('This script expects a single FASTA entry in the input file. Please try again.')
        
        fasta_lines = [line.replace('\n', '') for line in text if '>' not in line]
        fasta = ''.join(fasta_lines)
    
    if (prot_len := len(fasta)) > 400:
        raise Exception(f'The input protein is {prot_len} AA long.\nESMFold API query only allows proteins up to 400 AA long.\nTry using ColabFold instead.')
    
    # submit a new job via the API
    result = post('https://api.esmatlas.com/foldSequence/v1/pdb/',
                  data = fasta)
    
    with open(output_filepath, 'w+') as file:
        file.write(result._content.decode())

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
            
### NOTES
#ESMFold API example from website:
"""
curl -X POST --data "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL" https://api.esmatlas.com/foldSequence/v1/pdb/
"""
