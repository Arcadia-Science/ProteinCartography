#!/usr/bin/env python3
import argparse
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP
import numpy as np

MODES = ['pca', 'tsne', 'umap', 'pca_tsne', 'pca_umap']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = 'Path to input all-v-all similarity matrix.')
    parser.add_argument("-p", "--output-prefix", default = '', help = 'Prefix for resulting .tsv files.')
    parser.add_argument("-m", "--mode", default = 'pca', help = f'Mode of dimensionality reduction.\nValid arguments are {MODES}.')
    parser.add_argument("-r", "--random-state", default = '123456', help = 'Random state for umap and tsne modes.')
    args = parser.parse_args()
    
    return args

def calculate_PCA(pivot_file: str, n_components = 2, 
                  saveprefix = '', save = False, dimtype = 'pca', prep_step = False, **kwargs):
    pivoted_df = pd.read_csv(pivot_file, sep = '\t', index_col = 'protid')
    
    pca = PCA(n_components=n_components, **kwargs)
    pca_results = pca.fit_transform(pivoted_df)
    
    pca_results_df = pd.DataFrame(pca_results, columns = [f'PC{i}' for i in range(pca_results.ndim)], index = pivoted_df.index)
    
    if saveprefix != '':
        savefile = '_'.join([saveprefix, dimtype + '.tsv'])
    else:        
        savefile = pivot_file.replace('.tsv', '_' + dimtype + '.tsv')
    
    if save:
        pca_results_df.to_csv(savefile, sep = '\t')
        
    if prep_step:
        return savefile
    else:
        return pca_results_df
    
def calculate_TSNE(pivot_file: str, random_state: int,
                   n_components = 2, perplexity = 50, n_iter = 2000, 
                   saveprefix = '', save = False, dimtype = 'tsne', **kwargs):
    
    pivoted_df = pd.read_csv(pivot_file, sep = '\t', index_col = 'protid')
    
    if perplexity > len(pivoted_df):
        perplexity_check = int(1 + np.round(len(pivoted_df) / 5))
    else:
        perplexity_check = perplexity
    
    tsne = TSNE(n_components=n_components, perplexity = perplexity_check, 
                n_iter = n_iter, random_state = random_state, **kwargs)
    tsne_results = tsne.fit_transform(pivoted_df)
    
    tsne_results_df = pd.DataFrame(tsne_results, columns = [f'tSNE{i + 1}' for i in range(tsne_results.ndim)], index = pivoted_df.index)
    
    if saveprefix != '':
        savefile = '_'.join([saveprefix, dimtype + '.tsv'])
    else:        
        savefile = pivot_file.replace('.tsv', '_' + dimtype + '.tsv')
    
    if save:
        tsne_results_df.to_csv(savefile, sep = '\t')
    
    return tsne_results_df

def calculate_UMAP(pivot_file: str, random_state: int,
                   n_components = 2, n_neighbors= 80, min_dist = 0.5, 
                   saveprefix = '', save = False, dimtype = 'umap', **kwargs):
    
    pivoted_df = pd.read_csv(pivot_file, sep = '\t', index_col = 'protid')
    
    if n_neighbors > len(pivoted_df):
        neighbors_check = int(1 + np.round(len(pivoted_df) / 5))
    else:
        neighbors_check = n_neighbors
    
    umap_fxn = UMAP(n_components=n_components, random_state = random_state,
                n_neighbors = neighbors_check, min_dist = min_dist, **kwargs)
    umap_results = umap_fxn.fit_transform(pivoted_df)
    
    umap_results_df = pd.DataFrame(umap_results, columns = [f'UMAP{i + 1}' for i in range(umap_results.ndim)], index = pivoted_df.index)
    
    if saveprefix != '':
        savefile = '_'.join([saveprefix, dimtype + '.tsv'])
    else:        
        savefile = pivot_file.replace('.tsv', '_' + dimtype + '.tsv')
    
    if save:
        umap_results_df.to_csv(savefile, sep = '\t')
    
    return umap_results_df

def main():
    args = parse_args()
    pivot_file = args.input
    saveprefix = args.output_prefix
    mode = args.mode.lower()
    
    try:
        random_state = int(args.random_state)
    except:
        random_state = 123456
    
    if mode not in MODES:
        raise Exception(f'{mode} provided is not valid.\nValid modes include {MODES}.')
        
    if mode == 'pca':
        calculate_PCA(pivot_file, saveprefix = saveprefix, save = True)
        
    elif mode == 'tsne':
        calculate_TSNE(pivot_file, random_state, saveprefix = saveprefix, save = True)
        
    elif mode == 'umap':
        calculate_UMAP(pivot_file, random_state, saveprefix = saveprefix, save = True)
        
    elif mode == 'pca_tsne':
        pca_results_file = calculate_PCA(pivot_file, saveprefix = 
                                         saveprefix, save = True, prep_step = True)
        calculate_TSNE(pca_results_file, random_state, saveprefix = saveprefix, save = True)
        
    elif mode == 'pca_umap':
        pca_results_file = calculate_PCA(pivot_file, saveprefix = 
                                         saveprefix, save = True, prep_step = True)
        calculate_UMAP(pca_results_file, random_state, saveprefix = saveprefix, save = True)

if __name__ == '__main__':
    main()
