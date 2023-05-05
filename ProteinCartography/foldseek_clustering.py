#!/usr/bin/env python3
import argparse
import pandas as pd
import subprocess
import os
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query-folder", required = True)
    parser.add_argument("-r", "--results-folder", required = True)
    args = parser.parse_args()
    
    return args

def run_foldseek_clustering(query_folder: str, results_folder: str,
                            temp_folder = '', 
                            distances_filename = 'all_by_all_tmscore.tsv',
                            cluster_filename = 'clusters.tsv',
                            cluster_mode = '0', similarity_type = '2'):
    
    query_Path = Path(query_folder)
    
    if temp_folder == '':
        temp_Path = query_Path / 'temp'
    else:
        temp_Path = Path(temp_folder)
    
    results_Path = Path(results_folder)
    
    for path in [temp_Path, results_Path]:
        if not os.path.exists(path):
            os.mkdir(path)
    
    db_prefix = temp_Path / 'temp_afdb'
    subprocess.run(['foldseek', 'createdb', query_Path, db_prefix])

    foldseek_out = temp_Path / 'all_by_all'
    foldseek_tmp = temp_Path / 'tmp'
    subprocess.run(['foldseek', 'search', db_prefix, db_prefix, foldseek_out, foldseek_tmp, '-a'])

    foldseek_tmscore = temp_Path / 'all_by_all_tmscore'
    subprocess.run(['foldseek', 'aln2tmscore', db_prefix, db_prefix, foldseek_out, foldseek_tmscore])

    foldseek_distancestsv = results_Path / distances_filename
    subprocess.run(['foldseek', 'createtsv', db_prefix, db_prefix, foldseek_tmscore, foldseek_distancestsv])
    
    foldseek_cluster = temp_Path / 'clu'
    subprocess.run(['foldseek', 'clust', db_prefix, foldseek_out, foldseek_cluster, 
                    '--cluster-mode', cluster_mode, 
                    '--similarity-type', similarity_type])
    
    foldseek_clustertsv = results_Path / cluster_filename
    subprocess.run(['foldseek', 'createtsv', db_prefix, db_prefix, foldseek_cluster, foldseek_clustertsv])
    
    return str(foldseek_distancestsv), str(foldseek_clustertsv)
    
def make_struclusters_file(foldseek_clustertsv: str, output_file: str):
    df = pd.read_csv(foldseek_clustertsv, sep = '\t', names = ['ClusterRep', 'protid'])
    df['ClusterRep'] = df['ClusterRep'].str.rstrip('.pdb')
    df['protid'] = df['protid'].str.rstrip('.pdb')
    
    df_merged = df.groupby('ClusterRep').agg({i: ('first' if i == 'ClusterRep' else lambda x: [i for i in x]) for i in df.columns}).reset_index(drop = True)
    max_chars = len(str(df_merged.index[-1]))
    SC_ids = 'SC' + pd.Series(df_merged.index).apply(lambda x: str(x).zfill(max_chars))
    df_merged.insert(0, 'StruCluster', SC_ids)
    df_merged.drop(columns = ['ClusterRep'], inplace = True)
        
    df_exploded = df_merged.explode('protid')
    df_exploded = df_exploded[['protid', 'StruCluster']]
    df_exploded.to_csv(output_file, sep = '\t', index = None)

def clean_foldseek_results(input_file: str, output_file: str):
    """Takes the FoldSeek output file and cleans it to only contain the similarity scores and to/from values.
    
    Args:
        input_file (str): input filepath
        output_file (str): output filepath
    """
    with open(output_file, 'w+') as file:
        process1 = subprocess.Popen(('sed', "s/ /\t/g", input_file), stdout=subprocess.PIPE)
        process2 = subprocess.call(('awk', '-F', "\t", """{print $1 "\t" $2 "\t" $3}"""), stdin=process1.stdout, stdout=file)
    
def pivot_results(input_file: str, output_file: str):
    """Takes a df containing cleaned foldseek results and creates a similarity matrix.
    """
    import pandas as pd
    import os
    
    foldseek_df = pd.read_csv(input_file, sep = '\t', names = ['protid', 'target', 'tmscore'])
    pivoted_table = pd.pivot(foldseek_df, index = 'protid', columns = 'target', values = 'tmscore').fillna(0)
    
    pivoted_table.to_csv(output_file, sep = '\t')

def main():
    args = parse_args()
    query_folder = args.query_folder
    results_folder = args.results_folder
    
    distancestsv, clusterstsv = run_foldseek_clustering(query_folder, results_folder)
    
    cleanedtsv = distancestsv.replace('.tsv', '_cleaned.tsv')
    clean_foldseek_results(distancestsv, cleanedtsv)
    
    pivotedtsv = distancestsv.replace('.tsv', '_pivoted.tsv')
    pivot_results(cleanedtsv, pivotedtsv)
    
    featurestsv = clusterstsv.replace('.tsv', '_features.tsv')
    make_struclusters_file(clusterstsv, featurestsv)
    
if __name__ == '__main__':
    main()
    