#!/usr/bin/env python3
import argparse
import pandas as pd
import subprocess

BLAST_DEFAULTS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sacc', 'saccver', 'sgi', 'staxids', 'scomnames']
BLAST_DEFAULT_STRING = ' '.join(['6'] + BLAST_DEFAULTS)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = 'path of input peptide FASTA file.')
    parser.add_argument("-b", "--blastresults_output", required = True, help = 'path of destination blastresults.tsv file.')
    parser.add_argument("-o", "--output", required = True, help = 'path of destination blasthits.txt file, where blast hit accessions will be printed, one per line.')
    parser.add_argument("-M", "--max-target-seqs", default = '50000', help = 'maximum number of return hits. Default 50000')
    parser.add_argument("-B", "--blast-format-string", default = BLAST_DEFAULT_STRING, help = f'BLAST query format string.\n Default {BLAST_DEFAULT_STRING}')
    args = parser.parse_args()
    
    return args

def run_blast(input_file: str, output_file: str, 
              max_target_seqs = '50000', blast_format_string = BLAST_DEFAULT_STRING):
    '''
    Runs blastp using an input FASTA file and saves results to output.
    
    Args:
        input_file (str): path of input peptide FASTA file.
        output_file (str): path of destination blastresults.tsv file.
        max_target_seqs (str): maximum number of hits to return. Default '50000'.
        blast_format_string (str): passed to blastp '-outfmt'.
    '''
    subprocess.run([
        'blastp', 
        '-db', 'nr', 
        '-query', input_file, 
        '-out', output_file, 
        '-remote', 
        '-max_target_seqs', max_target_seqs, 
        '-outfmt', blast_format_string
    ])

# take an input blastresults file and create a .txt file from that
def extract_blasthits(input_file: str, output_file: str, 
                      names = BLAST_DEFAULTS):
    '''
    Takes an input blastresults.tsv file, reads the accessions, and prints unique hits to a .txt file, one per line.
    
    Args:
        input_file (str): path of input blastresults.tsv file.
        output_file (str): path of destination blasthits.txt file.
        names (str): names of columns of blast results.
    '''
    df = pd.read_csv(input_file, sep = '\t', names = names)

    hits = df['sacc'].unique()

    with open(output_file, 'w+') as f:
        f.writelines(hit + '\n' for hit in hits)
    
# run this if called from the interpreter
def main():
    args = parse_args()
    input_file = args.input
    blastresults_output = args.blastresults_output
    output_file = args.output
    
    max_target_seqs = args.max_target_seqs
    blast_format_string = args.blast_format_string
    
    if blast_format_string == BLAST_DEFAULT_STRING:
        blast_format_list = BLAST_DEFAULTS
    else:
        blast_format_list = [i for i in blast_format_string.split(' ') if i != '6']
    
    run_blast(input_file, blastresults_output, 
              max_target_seqs = max_target_seqs, 
              blast_format_string = blast_format_string)
    extract_blasthits(blastresults_output, output_file, 
                      names = blast_format_list)

# check if called from interpreter
if __name__ == '__main__':
    main()