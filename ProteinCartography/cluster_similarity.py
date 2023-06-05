#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import plotly.express as px
import arcadia_pycolor as apc

__all__ = ['calculate_group_similarity', 'plot_group_similarity']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix-file", required = True, help = 'Path to input all-v-all similarity matrix.')
    parser.add_argument("-f", "--features-file", required = True, help = 'Path to features file for grouping.')
    parser.add_argument("-c", "--features-column", required = True, help = 'Column to aggregate groups on.')
    parser.add_argument("-t", "--output-tsv", default = '', help = 'Path to output TSV file.')
    parser.add_argument("-h", "--output-html", default = '', help = 'Path to output HTML file.')
    args = parser.parse_args()
    
    return args

def calculate_group_similarity(matrix_file: str, features_file: str, features_column: str, output_file = None):
    '''
    '''
    # load in each dataframe
    pivot_df = pd.read_csv(matrix_file, sep = '\t')
    features_df = pd.read_csv(features_file, sep = '\t')
    
    # merge the two dataframes on protid
    pivot_agg = pivot_df.merge(features_df, on = 'protid')
    
    # group entries along one axis by feature groups, applying mean
    pivot_agg = pivot_agg.groupby(features_column).agg({i: list if i == 'protid' or i == features_column else np.mean for i in pivot_agg.columns})
    
    # tidy up dataframe
    pivot_agg.reset_index(drop = True, inplace = True)
    pivot_agg.drop(columns = ['protid', features_column], inplace = True)
    
    # transpose dataframe
    pivot_t = pivot_agg.transpose()
    pivot_t = pivot_t.reset_index().rename(columns = {'index': 'protid'})
    
    # merge features dataframe again and groupby feature
    pivot_t_agg = pivot_t.merge(features_df, on = 'protid')
    pivot_t_agg = pivot_t_agg.groupby(features_column).agg({i: list if i == 'protid' or i == features_column else np.mean for i in pivot_t_agg.columns})
    
    # clean up and reset groupings
    pivot_t_agg.drop(columns = ['protid', features_column], inplace = True)
    pivot_t_agg.columns = pivot_t_agg.index
    
    if output_file is not None:
        pivot_t_agg.to_csv(output_file, sep = '\t')
    
    return pivot_t_agg
    
def plot_group_similarity(group_similarity_file: str, plot_width = 700, plot_height = 700, output_file = None, show = False):
    arcadia_viridis = [(0, "#341E60"), 
                   (0.49, apc.arcadia_all["arcadia:aegean"]), 
                   (0.75, apc.arcadia_all["arcadia:lime"]),
                   (1, "yellow")
                  ]
    
    sim_df = pd.read_csv(group_similarity_file, index_col = 0, sep = '\t')

    fig = px.imshow(sim_df, color_continuous_scale = arcadia_viridis)
    fig.update_layout(width = plot_width, height = plot_height, 
                          coloraxis_colorbar=dict(
                              title="cluster similarity"
                          )
                     )
    fig.update_xaxes(side = "top")

    if output_file is not None:
        fig.write_html(output_file)
    
    if show:
        fig.show()
        
    return fig

def main():
    args = parse_args()
    matrix = args.matrix_file
    features = args.features_file
    column = args.features_column
    output_tsv = args.output_tsv
    output_html = args.output_html
    
    calculate_group_similarity(matrix, features, column, output_file = output_tsv)
    plot_group_similarity(output_tsv, output_file = output_html)
    
if __name__ == '__main__':
    main()
