#!/usr/bin/env python
import sys
import os
import argparse
from requests import get, post
from time import sleep

### NOTES
#FoldSeek API example from website:
"""
curl -X POST -F q=@PATH_TO_FILE -F 'mode=3diaa' \
-F 'database[]=afdb50' -F 'database[]=afdb-swissprot' -F 'database[]=afdb-proteome' \
-F 'database[]=mgnify_esm30' -F 'database[]=pdb100' -F 'database[]=gmgcl_id' \
https://search.foldseek.com/api/ticket
"""

# only import these functions when using import *
__all__ = ["foldseek_apiquery"]

# Possible align mode options from API
SET_MODES = ['3diaa', 'tmalign']

# Possible databases options from API
SET_DATABASES = ['afdb50', 'afdb-swissprot', 'afdb-proteome', 'mgnify_esm30', 'pdb100', 'gmgcl_id']
DEFAULT_DATABASES = ['afdb50', 'afdb-swissprot', 'afdb-proteome']

# parse command line arguments
def parse_args():
    # Set command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = 'Name of input file. Must end in ".pdb"')
    parser.add_argument("-o", "--output", required = True, 
                        help = 'Name of output file. Expects a ".tar.gz" suffix and will append if missing.')
    parser.add_argument("-m", "--mode", nargs = '?', default = '3diaa', 
                        help = ' | '.join([f"'{mode}'" for mode in SET_MODES]))
    parser.add_argument("-d", "--database", nargs = '+', default = DEFAULT_DATABASES, 
                        help = "'all' or any of " + ' | '.join([f"'{db}'" for db in SET_DATABASES]))
    args = parser.parse_args()
    return args

def foldseek_apiquery(input_file: str, output_file: str, mode: str, database: list):
    '''
    Queries the Foldseek web API with a PDB file and retrieves the results.
    
    Args:
        input_file (str): path to the query PDB file.
        output_file (str): path to a compressed '.tar.gz' results file. If suffix is missing, adds it.
        mode (str): whether to run in '3diaa' or 'tmalign' mode.
        database (list): list of run databases. 
            Valid databases include 'afdb50', 'afdb-swissprot', 'afdb-proteome', 'mgnify_esm30', 'pdb100', and 'gmgcl_id'.
    '''
    # Check to make sure input file has '.pdb' suffix
    if '.pdb' not in input_file:
        sys.exit('Input expects a .pdb file.')

    # Makes sure that the input file exists
    if not os.path.exists(input_file):
        sys.exit(f'File {input_file} not found.')

    # Append '.tar.gz' to file if it's not included
    if not output_file.endswith('.tar.gz'):
        output_file = output_file + '.tar.gz'

    # Checks for correct mode input
    if mode not in SET_MODES:
        sys.exit(f'Mode {mode} is not available. Accepted modes are {SET_MODES}.')

    ### Parse databases
    # Collector for user input databases
    query_databases = []

    # If all, use all the set databases
    if 'all' in database:
        query_databases = SET_DATABASES
        print(f'Querying all of the following databases: {query_databases}')
        
    # Otherwise, check to make sure each input database is valid
    else:
        for db in database:
            # Notify user that input database is not valid
            if db not in SET_DATABASES:
                print(f'{db} is not a valid option. ignoring.')
            else:
                query_databases.append(db)

    # Check to make sure at least one valid database is provided
    if len(query_databases) == 0:
        sys.exit(f'No valid databases provided. Valid databases include {SET_DATABASES}.')

    # Collector for PDB information for requests.post()
    pdb = ''

    # Open input file and collect text as string
    with open(input_file, 'r') as file:
        text = file.readlines()
        pdb = ''.join(text)
    
    ### Code below is mostly based on:
    ### <https://github.com/soedinglab/MMseqs2-App/blob/master/docs/api_example.py>
    # submit a new job via the API
    post_successful = False
    
    ticket = post('https://search.foldseek.com/api/ticket', {
                'q' : pdb,
                'database[]' : query_databases,
                'mode' : mode
            }).json()
        
    # check to see if the ticket failed to be posted
    # tickets can fail to be posted because of ratelimits
    try:
        testid = ticket['id']
    except KeyError:
        if 'status' in ticket.keys() and 'reason' in ticket.keys():
            print('===============')
            print(ticket['status'])
            print(ticket['reason'])
            print('===============')
        sys.exit('Foldseek may be rate-limiting your requests. Try again later.')

    # poll until the job was successful or failed
    repeat = True
    tries = 0
    limit = 10
    sleep_time = 30
    while repeat and tries < limit:
        status = get('https://search.foldseek.com/api/ticket/' + ticket['id']).json()
        if status['status'] == "ERROR":
        # handle error
            sys.exit('The ticket returned with status ERROR.')

        # wait a short time between poll requests
        sleep(sleep_time)
        tries += 1
        repeat = status['status'] != "COMPLETE"
    
    if tries == 10:
        sys.exit(f'The ticket failed to complete after {tries * sleep_time} seconds.')

    # download blast compatible result archive
    download = get('https://search.foldseek.com/api/result/download/' + ticket['id'], stream=True)
    with open(output_file, 'wb') as fd:
        for chunk in download.iter_content(chunk_size=128):
            fd.write(chunk)

# run this if called from the interpreter
def main():
    # parse args
    args = parse_args()
    
    input_file = args.input
    output_file = args.output
    mode = args.mode
    database = args.database
    
    foldseek_apiquery(input_file, output_file, mode, database)

# check if called from interpreter
if __name__ == '__main__':
    main()
