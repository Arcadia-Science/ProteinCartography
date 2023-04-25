#!/usr/bin/env python3
import sys
import os
import argparse
from requests import get, post
from time import sleep

# Possible align mode options from API
set_modes = ['3diaa', 'tmalign']

# Possible databases options from API
set_databases = ['afdb50', 'afdb-swissprot', 'afdb-proteome', 'mgnify_esm30', 'pdb100', 'gmgcl_id']

# Set command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required = True, help = 'Name of input file. Must end in ".pdb"')
parser.add_argument("-o", "--output", required = True, help = 'Name of output file. Expects a ".tar.gz" suffix and will append if missing.')
parser.add_argument("-m", "--mode", nargs = '?', default = '3diaa', help = ' | '.join([f"'{mode}'" for mode in set_modes]))
parser.add_argument("-d", "--database", nargs = '+', default = ['afdb50', 'afdb-swissprot', 'afdb-proteome'], help = "'all' or any of " + ' | '.join([f"'{db}'" for db in set_databases]))
args = parser.parse_args()

# Check to make sure input file has '.pdb' suffix
if '.pdb' not in (input_file := args.input):
    sys.exit('Input expects a .pdb file.')

# Makes sure that the input file exists
if not os.path.exists(input_file):
    sys.exit(f'File {input_file} not found.')

# Append '.tar.gz' to file if it's not included
if '.tar.gz' not in (output_file := args.output):
    print('appending ".tar.gz" to output name')
    
    output_file = output_file + '.tar.gz'

# Checks for correct mode input
if (mode := args.mode) not in set_modes:
    sys.exit(f'Mode {mode} is not available. Accepted modes are {set_modes}.')

### Parse databases
# Collector for user input databases
query_databases = []

# If all, use all the set databases
if 'all' in args.database:
    query_databases = set_databases
    print(f'Querying all of the following databases: {query_databases}')
# Otherwise, check to make sure each input database is valid
else:
    for db in args.database:
        # Notify user that input database is not valid
        if db not in set_databases:
            print(f'{db} is not a valid option. ignoring.')
        else:
            query_databases.append(db)

# Check to make sure at least one valid database is provided
if len(query_databases) == 0:
    sys.exit(f'No valid databases provided. Valid databases include {set_databases}.')

# Collector for PDB information for requests.post()
pdb = ''

# Open input file and collect text as string
with open(input_file, 'r') as file:
    text = file.readlines()
    pdb = ''.join(text)
    
### Code below is mostly based on:
### 
# submit a new job via the API
ticket = post('https://search.foldseek.com/api/ticket', {
            'q' : pdb,
            'database[]' : query_databases,
            'mode' : mode
        }).json()

# poll until the job was successful or failed
repeat = True
while repeat:
    status = get('https://search.foldseek.com/api/ticket/' + ticket['id']).json()
    if status['status'] == "ERROR":
        # handle error
        sys.exit(0)

    # wait a short time between poll requests
    sleep(1)
    repeat = status['status'] != "COMPLETE"

# download blast compatible result archive
download = get('https://search.foldseek.com/api/result/download/' + ticket['id'], stream=True)
with open(output_file, 'wb') as fd:
    for chunk in download.iter_content(chunk_size=128):
        fd.write(chunk)

        
### NOTES
#FoldSeek API example from website:
"""
curl -X POST -F q=@PATH_TO_FILE -F 'mode=3diaa' \
-F 'database[]=afdb50' -F 'database[]=afdb-swissprot' -F 'database[]=afdb-proteome' \
-F 'database[]=mgnify_esm30' -F 'database[]=pdb100' -F 'database[]=gmgcl_id' \
https://search.foldseek.com/api/ticket
"""
