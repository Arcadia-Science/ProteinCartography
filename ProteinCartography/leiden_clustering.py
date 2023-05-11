import argparse
import scanpy as sc
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = 'Input file path of a similarity matrix.')
    parser.add_argument("-o", "--output", required = True, help = 'Output path to a file, usually leiden_features.tsv')
    parser.add_argument("-n", "--neighbors", default = '10', help = 'Number of n_neighbors to pass to sc.pp.neighbors().')
    parser.add_argument("-c", "--components", default = '30', help = 'Number of n_pcs to pass to sc.pp.neighbors().')
    args = parser.parse_args()
    return args

def scanpy_leiden_cluster(input_file: str, savefile = '', n_neighbors = 10, n_pcs = 30):
    adata = sc.read_csv(input_file, delimiter = '\t')
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, n_pcs = n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    membership = pd.DataFrame(adata.obs['leiden']).reset_index()
    membership.rename(columns = {'index': 'protid', 'leiden': 'LeidenCluster'}, inplace = True)
    max_chars = len(str(membership['LeidenCluster'].astype(int).max()))
    membership['LeidenCluster'] = 'LC' + membership['LeidenCluster'].apply(lambda x: str(x).zfill(max_chars)).astype(str)
    
    if savefile != '':
        membership.to_csv(savefile, sep = '\t', index = None)
    
    return membership

def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    neighbors = int(args.neighbors)
    pcs = int(args.components)
    
    scanpy_leiden_cluster(input_file = input_file, savefile = output_file, n_neighbors = neighbors, n_pcs = pcs)
    
if __name__ == '__main__':
    main()
