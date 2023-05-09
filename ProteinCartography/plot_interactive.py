import argparse
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import arcadia_pycolor as apc
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dimensions", required = True)
    parser.add_argument("-f", "--features", required = True)
    parser.add_argument("-o", "--output", required = True)
    args = parser.parse_args()
    
    return args

def apply_coordinates(dimensions_file: str, features_file: str, saveprefix = '', save = False, prep_step = False):
    reduced_dim_df = pd.read_csv(dimensions_file, sep = '\t')
    dimtype = ''.join([i for i in reduced_dim_df.columns[1] if not i.isdigit()])
    agg_features_df = pd.read_csv(features_file, sep = '\t')
    
    if '.tsv' not in features_file:
        raise Exception(f'{features_file} does not end in ".tsv" as expected.')

    plot_features_df = reduced_dim_df.merge(agg_features_df, on = 'protid')
    
    if saveprefix != '':
        savefile = '_'.join([saveprefix, dimtype + '.tsv'])
    else:        
        savefile = features_file.replace('.tsv', '_' + dimtype + '.tsv')
    
    if save:
        plot_features_df.to_csv(savefile, sep = '\t', index = None)
    
    if prep_step:
        return savefile
    
    return plot_features_df

def assign_taxon(taxon_list: list, rank_list: list, hierarchical = False, sep = ','):
    # check for every element of the rank list whether it's in the taxon list
    output = [item for item in rank_list if item in taxon_list]
    
    # if output is empty, fill it in with something
    if output == []:
        output.append('Other')
    
    # if hierarchical, return the lowest rank item
    if hierarchical:
        return output[0]
    # else, return a joined string
    else:
        return sep.join(output)

def plot_interactive(coordinates_file: str, plotting_rules: dict,
                    marker_size = 5, marker_opacity = 0.8, output_file = '', show = False):
    
    # make a copy of the starting data
    df = pd.read_csv(coordinates_file, sep = '\t')
    
    # generate a hover template text string
    hovertemplate_generator = ["<b>%{customdata[0]}</b></br>––––––––––––"]
    # pass in the custom data columns
    custom_data = ['protid'] + list(plotting_rules.keys())
    
    # iterate through plotting rules, applying relevant rules
    for col in plotting_rules:
        # if the plotting rule 'fillna' is present, fills NAs with that value
        if 'fillna' in plotting_rules[col].keys():
            df[col] = df[col].fillna(plotting_rules[col]['fillna'])
            
        # if the plotting rule 'apply' is present, applies that function
        if 'apply' in plotting_rules[col].keys():
            df[col] = df[col].apply(plotting_rules[col]['apply'])
        
        # if you want to skip an element from being added to hover text
        if 'skip_hover' in plotting_rules[col].keys():
            continue
        
        # sets up hover text based on 'textlabel' attributes
        if 'textlabel' in plotting_rules[col].keys():
            hovertext_label = plotting_rules[col]['textlabel']
        else:
            hovertext_label = col
            
        # generates hoverlabel custom text string for that column
        # This doesn't work properly if the column is NA for that value.
        col_index = custom_data.index(col)
        
        hovertemplate_item = "<b>" + hovertext_label + "</b>: %{customdata[" + str(col_index) + "]}"
        hovertemplate_generator.append(hovertemplate_item)
        
    # generates a full hovertemplate string from hovertemplate_generator list
    hovertemplate = "<br>".join(hovertemplate_generator)
    
    # gets the first two PCs or tSNE columns or UMAP columns
    dim1 = df.columns[1]
    dim2 = df.columns[2]
    
    # Collector dictionary for making plots
    plots = {}
    
    # Iterate through the plotting rules for each plotted datatype
    # We make each plot for each datatype individually;
    # Then, we transfer the points from the existing plots to a new, single plot
    # This way we can get the toggle system working
    for col in plotting_rules.keys():
        # Plotting rules for categorical plots
        if plotting_rules[col]['type'] == 'categorical':
            
            color_keys = np.sort(df[col].unique())
            
            # if a color dict is specified, use that color dict
            if 'color_dict' in plotting_rules[col].keys():
                
                colors_dict = plotting_rules[col]['color_dict']
                
                # make sure there's at least one entry for every possible value for that column
                # in the future, should automatically add extra colors.
                if not all([entry in colors_dict.keys() for entry in color_keys]):
                    raise Exception('color_dict is missing some entries.')
            
            # if a color dict is not provided, but a color order is, use that to make a color dict
            elif 'color_order' in plotting_rules[col].keys():
                color_order = plotting_rules[col]['color_order']
                
                # make sure there's enough colors to make a dict with
                # in the future, should automatically add extra colors
                if len(color_order) < len(color_keys):
                    raise Exception(f'color_order expects an equal or greater number of colors for unique values.\nProvided {len(color_order)} colors for {len(color_keys)} values.')
                colors_dict = dict(zip(color_keys, color_order))
            
            else:
                colors_dict = dict(zip(color_keys))
            
            # generate a plot using the above parameters
            plots[col] = px.scatter(df, dim1, dim2, color = col, hover_name = 'protid', 
                                   category_orders = {col: color_keys},
                                   color_discrete_map = colors_dict,
                                   custom_data = custom_data)
            # add the hovertemplate text and other aspects to the traces for that column
            plots[col].update_traces(marker = dict(size = marker_size, opacity = marker_opacity),
                                    hovertemplate = hovertemplate)
        
        # Plotting rules for continuous plots
        elif plotting_rules[col]['type'] == 'continuous':
            
            # Color scales seem to be broken for some reason; they always show as Magma.
            # When plotting a single plot, this works fine.
            # However, moving the points to a new plot breaks the color scheme.
            # Try to find a workaround in the future.
            if 'color_scale' in plotting_rules[col].keys():
                color_scale = plotting_rules[col]['color_scale']
            else:
                color_scale = 'viridis'
            
            # generate a plot using the above parameters
            plots[col] = px.scatter(df, dim1, dim2, color = col, hover_name = 'protid',
                                   color_continuous_scale = color_scale,
                                   custom_data = custom_data)
            plots[col].update_traces(marker = dict(size = marker_size, opacity = marker_opacity),
                                    hovertemplate = hovertemplate)
        
        elif plotting_rules[col]['type'] == 'taxonomic':
            
            if 'taxon_order' in plotting_rules[col].keys():
                taxon_order = plotting_rules[col]['taxon_order']
            else:
                raise Exception('Please provide a "taxon_order" list in the plotting rules.')
            
            color_keys = taxon_order
            
            # if a color dict is not provided, but a color order is, use that to make a color dict
            if 'color_order' in plotting_rules[col].keys():
                color_order = plotting_rules[col]['color_order']
                
                # make sure there's enough colors to make a dict with
                # in the future, should automatically add extra colors
                if len(color_order) < len(color_keys):
                    raise Exception(f'color_order expects an equal or greater number of colors for unique values.\nProvided {len(color_order)} colors for {len(color_keys)} values.')
                colors_dict = dict(zip(color_keys, color_order))
                colors_dict['Other'] = '#eaeaea'
            
            col_taxonomic = col + '_taxonomic'
            df[col_taxonomic] = df[col].apply(lambda x: assign_taxon(x, taxon_order, hierarchical = True))
            
            plots[col] = px.scatter(df, dim1, dim2, color = col_taxonomic, hover_name = 'protid',
                                   color_discrete_map = colors_dict,
                                   custom_data = custom_data)
            plots[col].update_traces(marker = dict(size = marker_size, opacity = marker_opacity),
                                    hovertemplate = hovertemplate)
            
    scatter_counter = {}
    fig_order = []
    
    fig = go.Figure()
    for j, (col, plot) in enumerate(plots.items()):
        if j == 0:
            vis = True
        else:
            vis = False
            
        fig_order.append(col)
        
        scatter_counter[col] = 0        
        for scatter in plot.data:
            fig.add_trace(go.Scatter(scatter, visible = vis))
            scatter_counter[col] += 1
    
    def visibility_list(col):
        entries = []
        
        for fig in fig_order:
            if fig == col:
                entries.extend([True] * scatter_counter[fig])
            else:
                entries.extend([False] * scatter_counter[fig])
        
        return entries
    
    buttons = []
    for k, (col, plot) in enumerate(plots.items()):
        if 'textlabel' in plotting_rules[col].keys():
            button_label = plotting_rules[col]['textlabel']
        else:
            button_label = col
        
        button_item = dict(
                        args = ["visible", visibility_list(col)],
                        label = button_label,
                        method = "restyle"
                    )
        buttons.append(button_item)

    fig.update_layout(
        updatemenus=[
            dict(
                buttons = list(buttons),
                showactive = True,
                x = 0.1,
                xanchor = "left",
                y = 1.1,
                yanchor = "top"
            ),
        ]
    )

    fig.update_layout(
        width = 600, height = 600,
        annotations=[
            dict(text="color", x=0.01, xref="paper", y=1.08, yref="paper",
             align="left", showarrow=False),
        ],
        plot_bgcolor= '#fafafa'
    )
    fig.update_layout(xaxis = dict(showticklabels = False), yaxis = dict(showticklabels = False))
    fig.update_yaxes(scaleanchor = "x", scaleratio = 1)
    fig.update_layout(legend=dict(orientation = "h", yanchor="top", y = 0, xanchor="left", x = 0))
    fig.update_layout(coloraxis_colorbar=dict(yanchor="top", y=0, x=0.5,
                                          ticks="outside", orientation = 'h'))
    if output_file != '':
        fig.write_html(output_file)
    if show:
        fig.show()

def main():
    args = parse_args()
    dimensions_file = args.dimensions
    features_file = args.features
    output_file = args.output

    annotationScore_colors = [
        apc.arcadia_all['arcadia:crow'], 
        apc.arcadia_all['arcadia:aster'],
        apc.arcadia_all['arcadia:aegean'],
        apc.arcadia_all['arcadia:seaweed'],
        apc.arcadia_all['arcadia:lime'],
        apc.arcadia_all['arcadia:canary']
    ]

    annotationScore_color_dict = dict(zip([str(i) for i in range(6)], annotationScore_colors))
    taxon_color_dict = {
        'Mammalia': apc.arcadia_all['arcadia:oat'],
        'Vertebrata': apc.arcadia_all['arcadia:canary'],
        'Arthropoda': apc.arcadia_all['arcadia:seaweed'],
        'Ecdysozoa': apc.arcadia_all['arcadia:mint'],
        'Lophotrocozoa': apc.arcadia_all['arcadia:aegean'],
        'Metazoa': apc.arcadia_all['arcadia:amber'],
        'Fungi': apc.arcadia_all['arcadia:periwinkle'], 
        'Eukaryota': apc.arcadia_all['arcadia:aster'], 
        'Bacteria': apc.arcadia_all['arcadia:slate'], 
        'Archaea': apc.arcadia_all['arcadia:dragon']
    }

    plotting_rules = {
        'organism.scientificName': {
            'type': 'hovertext',
            'fillna': '',
            'textlabel': 'Species'
        },
        'organism.commonName': {
            'type': 'hovertext',
            'fillna': '',
            'textlabel': 'Common name'
        },
        'StruCluster': {
            'type': 'categorical',
            'fillna': 'None',
            'color_order': apc.arcadia_All_ordered.values(),
            'textlabel': 'Structural Cluster'
        },
        'annotationScore': {
            'type': 'categorical',
            'fillna': 0,
            'apply': lambda x: str(int(x)),
            'color_dict': annotationScore_color_dict,
            'textlabel': 'Annotation Score'
        },
        'organism.lineage': {
            'type': 'taxonomic',
            'fillna': '[]',
            'apply': lambda x: eval(x),
            'taxon_order': taxon_color_dict.keys(),
            'color_order': taxon_color_dict.values(),
            'textlabel': 'Broad Taxon',
            'skip_hover': True
        },
        'sequence.length': {
            'type': 'continuous',
            'fillna': 0,
            'textlabel': 'Length'
        }
    }
    
    coordinates_file = apply_coordinates(dimensions_file, features_file, save = True, prep_step = True)
    plot_interactive(coordinates_file, plotting_rules, output_file = output_file)
    
if __name__ == '__main__':
    main()
