#!/usr/bin/env python
import argparse
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP
import numpy as np

# only import these functions when using import *
__all__ = ["calculate_PCA", "calculate_TSNE", "calculate_UMAP"]

MODES = ['pca', 'tsne', 'umap', 'pca_tsne', 'pca_umap']

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = 'Path to input all-v-all similarity matrix.')
    parser.add_argument("-p", "--output-prefix", default = '', help = 'Prefix for resulting .tsv files.')
    parser.add_argument("-m", "--mode", default = 'pca', help = f'Mode of dimensionality reduction.\nValid arguments are {MODES}.')
    parser.add_argument("-r", "--random-state", default = '123456', help = 'Random state for umap and tsne modes.')
    args = parser.parse_args()
    
    return args

def calculate_PCA(pivot_file: str, n_components = 2, 
                  save = False, saveprefix = None, dimtype = 'pca', prep_step = False, **kwargs):
    '''
    Calculates the n-component PCA of an all-v-all similarity matrix.
    The matrix should be a square, where rows and columns are identical and each cell is a similarity score.
    
    Args:
        pivot_file (str): path to a matrix of values.
        n_components (int): number of components to calculate. Default 2.
        save (bool): whether or not to save the file.
        saveprefix (str): prefix of file to save to.
        dimtype (str): defaults to 'pca'. included in save file output name.
        prep_step (bool): if set to True, returns the path of the saved file.
        **kwargs are passed to `sklearn.decomposition.PCA`.
    Returns:
        a pandas.DataFrame containing the PCA results, or a path to the saved file.
    '''
    # Read input file
    pivoted_df = pd.read_csv(pivot_file, sep = '\t', index_col = 'protid')
    
    # Initialize and run PCA
    pca = PCA(n_components=n_components, **kwargs)
    pca_results = pca.fit_transform(pivoted_df)
    
    # Read PCA results data
    pca_results_df = pd.DataFrame(pca_results, columns = [f'PC{i}' for i in range(pca_results.shape[1])], index = pivoted_df.index)
    
    # Generate savefile name
    if saveprefix is not None:
        savefile = '_'.join([saveprefix, dimtype + '.tsv'])
    else:        
        savefile = pivot_file.replace('.tsv', '_' + dimtype + '.tsv')
    
    # Save if needed
    if save:
        pca_results_df.to_csv(savefile, sep = '\t')
        
    # Return either the path to the destination file or the results dataframe
    if prep_step:
        return savefile
    else:
        return pca_results_df
    
def calculate_TSNE(pivot_file: str, random_state: int,
                   n_components = 2, perplexity = 50, n_iter = 2000, 
                   save = False, saveprefix = None, dimtype = 'tsne', **kwargs):
    '''
    Calculates the TSNE of an all-v-all similarity matrix given a random state, perplexity, and number of iterations.
    The matrix should be a square, where rows and columns are identical and each cell is a similarity score.
    
    Args:
        pivot_file (str): path to a matrix of values.
        random_state (int): random state used for initializing TSNE.
        n_components (int): number of components to return. Default 2.
        perplexity (int): `sklearn.manifold.TSNE` perplexity.
        n_iter (int): iteractions to run with TSNE.
        save (bool): whether or not to save the file.
        saveprefix (str): prefix of file to save to.
        dimtype (str): defaults to 'tsne'. included in save file output name.
        **kwargs are passed to `sklearn.manifold.TSNE`.
    Returns:
        a pandas.DataFrame containing the TSNE results, or a path to the saved file.
    '''
    # Read input file
    pivoted_df = pd.read_csv(pivot_file, sep = '\t', index_col = 'protid')
    
    # Check to make sure perplexity is lower than the total number of elements
    # If not, set perplexity to 1/5 of elements
    if perplexity > len(pivoted_df):
        perplexity_check = int(1 + np.round(len(pivoted_df) / 5))
    else:
        perplexity_check = perplexity
    
    # Intialize and run TSNE
    tsne = TSNE(n_components=n_components, perplexity = perplexity_check, 
                n_iter = n_iter, random_state = random_state, **kwargs)
    tsne_results = tsne.fit_transform(pivoted_df)
    
    # Read TSNE results data
    tsne_results_df = pd.DataFrame(tsne_results, columns = [f'tSNE{i + 1}' for i in range(tsne_results.ndim)], index = pivoted_df.index)
    
    # Generate savefile name
    if saveprefix is not None:
        savefile = '_'.join([saveprefix, dimtype + '.tsv'])
    else:        
        savefile = pivot_file.replace('.tsv', '_' + dimtype + '.tsv')
    
    # Save if needed
    if save:
        tsne_results_df.to_csv(savefile, sep = '\t')
    
    # Return results
    return tsne_results_df

def calculate_UMAP(pivot_file: str, random_state: int,
                   n_components = 2, n_neighbors= 80, min_dist = 0.5, 
                   save = False, saveprefix = None, dimtype = 'umap', **kwargs):
    '''
    Calculates the UMAP of an all-v-all similarity matrix given a random state, number of neighbors, and minimum distance.
    The matrix should be a square, where rows and columns are identical and each cell is a similarity score.
    
    Args:
        pivot_file (str): path to a matrix of values.
        random_state (int): random state used for initializing UMAP.
        n_components (int): number of components to return. Default 2.
        n_neighbors (int): number of neighbors.
        min_dist (int): minimum distance between neighbors.
        save (bool): whether or not to save the file.
        saveprefix (str): prefix of file to save to.
        dimtype (str): defaults to 'tsne'. included in save file output name.
        **kwargs are passed to `umap.UMAP`.
    Returns:
        a pandas.DataFrame containing the TSNE results, or a path to the saved file.
    '''
    # Read input file
    pivoted_df = pd.read_csv(pivot_file, sep = '\t', index_col = 'protid')
    
    # Check to make sure number of neighbors isn't greater than the whole dataset
    # If it is, set number of neighbors to 1/5 of data
    if n_neighbors > len(pivoted_df):
        neighbors_check = int(1 + np.round(len(pivoted_df) / 5))
    else:
        neighbors_check = n_neighbors
    
    # Initialize and run UMAP
    umap_fxn = UMAP(n_components=n_components, random_state = random_state,
                n_neighbors = neighbors_check, min_dist = min_dist, **kwargs)
    umap_results = umap_fxn.fit_transform(pivoted_df)
    
    # Read UMAP results data
    umap_results_df = pd.DataFrame(umap_results, columns = [f'UMAP{i + 1}' for i in range(umap_results.ndim)], index = pivoted_df.index)
    
    # Generate savefile name
    if saveprefix is not None:
        savefile = '_'.join([saveprefix, dimtype + '.tsv'])
    else:        
        savefile = pivot_file.replace('.tsv', '_' + dimtype + '.tsv')
    
    # Save if needed
    if save:
        umap_results_df.to_csv(savefile, sep = '\t')
    
    # Return results
    return umap_results_df

# run this if called from the interpreter
def main():
    args = parse_args()
    pivot_file = args.input
    saveprefix = args.output_prefix
    mode = args.mode.lower()
    
    # Set a random state based on user input
    # Or use default value
    try:
        random_state = int(args.random_state)
    except:
        random_state = 123456
    
    # Check whether modes are valid
    if mode not in MODES:
        raise Exception(f'{mode} provided is not valid.\nValid modes include {MODES}.')
        
    # Calculate vanilla PCA
    if mode == 'pca':
        calculate_PCA(pivot_file, save = True, saveprefix = saveprefix)
    
    # Calculate TSNE
    elif mode == 'tsne':
        calculate_TSNE(pivot_file, random_state, save = True, saveprefix = saveprefix)
    
    # Calcualte UMAP
    elif mode == 'umap':
        calculate_UMAP(pivot_file, random_state, save = True, saveprefix = saveprefix)
    
    # First run 30-component PCA and pass that space to TSNE
    elif mode == 'pca_tsne':
        # Generate temporary save prefix to avoid collisions if run simultaneously with other PCA calls
        saveprefix1 = pivot_file.replace('.tsv', 'temp1')
        pca_results_file = calculate_PCA(pivot_file, save = True, saveprefix = saveprefix1,
                                         n_components = 30, prep_step = True)
        
        saveprefix2 = pca_results_file.replace('temp1', '').replace('.tsv', '')
        calculate_TSNE(pca_results_file, random_state, save = True, saveprefix = saveprefix2)
    
    # First run 30-component PCA and pass that space to UMAP
    elif mode == 'pca_umap':
        # Generate temporary save prefix to avoid collisions if run simultaneously with other PCA calls
        saveprefix1 = pivot_file.replace('.tsv', 'temp2')
        pca_results_file = calculate_PCA(pivot_file, save = True, saveprefix = saveprefix1,
                                         n_components = 30, prep_step = True)
        
        saveprefix2 = pca_results_file.replace('temp2', '').replace('.tsv', '')
        calculate_UMAP(pca_results_file, random_state, save = True, saveprefix = saveprefix2)

# check if called from interpreter
if __name__ == '__main__':
    main()
