#!/usr/bin/env python
from bioservices import UniProt
import argparse
import os
import subprocess

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--accession", required = True, nargs = '+')
    parser.add_argument("-o", "--output", required = True)
    parser.add_argument("-f", "--format", nargs = '+', default = ['fasta', 'pdb'])
    args = parser.parse_args()
    return args

def fetch_fasta(accession: str, output_dir: str):
    u = UniProt()
    output_path = output_dir + accession + '.fasta'
        
    if not os.path.exists(output_path):
        res = u.retrieve(accession, frmt = 'fasta')
        with open(output_path, 'w+') as f:
            f.write(res)
    
def fetch_pdb(accession: str, output_dir: str):
    output_path = output_dir + accession + '.pdb'
    source = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'.format(accession)

    if not os.path.exists(output_path):
        subprocess.run(['curl' , '-JLo' , output_path, source], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
def main():
    args = parse_args()
    accessions = args.accession
    output_dir = args.output
    
    if output_dir[-1] != '/':
        output_dir = output_dir + '/'
    
    formats = args.format
    
    for accession in accessions:
        if 'fasta' in formats:
            fetch_fasta(accession, output_dir)
        if 'pdb' in formats:
            fetch_pdb(accession, output_dir)
    
if __name__ == '__main__':
    main()
