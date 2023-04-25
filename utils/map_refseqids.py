from bioservices import UniProt
import argparse
import pandas as pd
u = UniProt()

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input")
parser.add_argument("-o", "--output")
args = parser.parse_args()

input_file = args.input
output_file = args.output

query_dbs = ['EMBL-GenBank-DDBJ_CDS', 'RefSeq_Protein']

with open(input_file, 'r') as f:
    ids = f.read().splitlines()

dummy_df = pd.DataFrame()
    
for i, db in enumerate(query_dbs):
    results = u.mapping(db, "UniProtKB", query = ' '.join(ids))
    results_df = pd.json_normalize(results['results'])
    
    if i == 0:
        dummy_df = results_df
    else:
        dummy_df = pd.concat([dummy_df, results_df], axis = 0)

hits = dummy_df['to.primaryAccession'].unique()

with open(output_file, 'w+') as f:
    f.writelines(hit + '\n' for hit in hits)