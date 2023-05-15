import argparse
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import arcadia_pycolor as apc
import numpy as np
import matplotlib.colors as mc
import colorsys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dimensions", required = True)
    parser.add_argument("-f", "--features", required = True)
    parser.add_argument("-o", "--output", required = True)
    parser.add_argument("-t", "--dimensions-type", default = '')
    parser.add_argument("-k", "--keyids", nargs = "+", default = [])
    parser.add_argument("-x", "--taxon_focus", default = "euk")
    args = parser.parse_args()
    
    return args

def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    c2 = colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])
    return mc.to_hex(c2)

def apply_coordinates(dimensions_file: str, features_file: str, saveprefix = '', dimtype = '', save = False, prep_step = False):
    reduced_dim_df = pd.read_csv(dimensions_file, sep = '\t')
    if dimtype == '':
        outdimtype = ''.join([i for i in reduced_dim_df.columns[1] if not i.isdigit()])
    else:
        outdimtype = dimtype
    agg_features_df = pd.read_csv(features_file, sep = '\t')
    
    if '.tsv' not in features_file:
        raise Exception(f'{features_file} does not end in ".tsv" as expected.')

    plot_features_df = reduced_dim_df.merge(agg_features_df, on = 'protid')
    
    if saveprefix != '':
        savefile = '_'.join([saveprefix, outdimtype + '.tsv'])
    else:        
        savefile = features_file.replace('.tsv', '_' + outdimtype + '.tsv')
    
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
                    marker_size = 4, marker_opacity = 0.8, output_file = '', keyids = [], show = False):
    
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
                color_order = list(plotting_rules[col]['color_order'])
                
                # make sure there's enough colors to make a dict with
                # in the future, should automatically add extra colors
                if len(color_order) < len(color_keys):
                    num_cycles = int(np.round(len(color_keys) / len(color_order)))
                    
                    more_colors = []
                    steps = [0.7, 0.5, 0.3]
                    
                    if num_cycles > len(steps):
                        raise Exception(f'Can create up to {len(steps) * len(color_order)} colors to use.\nNeeded {len(color_keys)} colors.')
                    
                    for n in range(num_cycles - 1):
                        color_order_duplicated = [adjust_lightness(color) for color in color_order]
                        more_colors.extend(color_order_duplicated)
                    
                    color_order.extend(more_colors)
                    
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
            
            if len(colors_dict.keys()) > 20:
                plots[col].layout.update(showlegend = False)
        
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
        
        scatter_counter[col] = len(plot.data)
        if scatter_counter[col] > 15:
            showlegend = False
        else:
            showlegend = True
        
        for scatter in plot.data:
            fig.add_trace(go.Scattergl(scatter, visible = vis, showlegend = showlegend))
        
    if keyids != []:
        keypoints = df[df['protid'].isin(keyids)]
        
        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=keypoints[dim1],
                y=keypoints[dim2],
                marker=dict(
                    color='rgba(0,0,0,1)',
                    size=10,
                    symbol='star-diamond'
                ),
                showlegend=False,
                hoverinfo='skip',
                hoverlabel=None,
                name = 'keyids'
            )
        )
    
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
                        args = ["visible", visibility_list(col) + [True]],
                        label = button_label,
                        method = "restyle"
                    )
        buttons.append(button_item)
        
    dropdown_menu = dict(
                buttons = list(buttons),
                showactive = True,
                x = 0.1,
                xanchor = "left",
                y = 1.1,
                yanchor = "top"
    )
    
    if keyids != []:
        keyids_button = dict(
            buttons = [dict(
                args = ['visible', [i.visible if i.name != 'keyids' else True for i in fig.select_traces()]],
                args2 = ['visible', [i.visible if i.name != 'keyids' else False for i in fig.select_traces()]],
                label = 'Input Proteins'
            )],
            type = 'buttons',
            showactive = True,
            x = 0.5,
            xanchor = "left",
            y = 1.1,
            yanchor = "top" 
        )
        updatemenu = [dropdown_menu, keyids_button]
    else:
        updatemenu = [dropdown_menu]

    fig.update_layout(
        updatemenus = updatemenu
    )

    fig.update_layout(
        width = 700, height = 700,
        annotations=[
            dict(text="color", x=0.01, xref="paper", y=1.08, yref="paper",
             align="left", showarrow=False),
        ],
        plot_bgcolor= '#fafafa'
    )
    fig.update_layout(xaxis = dict(showticklabels = False), yaxis = dict(showticklabels = False))
    fig.update_yaxes(scaleanchor = "x", scaleratio = 1)
    fig.update_layout(legend=dict(orientation = "h", yanchor="top", y = 0, xanchor="left", x = 0, font = dict(size = 10)))
    fig.update_layout(coloraxis_colorbar=dict(yanchor="top", y=0, x=0.5,
                                          ticks="outside", orientation = 'h'))
    if output_file != '':
        fig.write_html(output_file)
    
    if show:
        fig.show()
    
    return fig

def main():
    args = parse_args()
    dimensions_file = args.dimensions
    features_file = args.features
    output_file = args.output
    dimensions_type = args.dimensions_type
    keyids = args.keyids
    taxon_focus = args.taxon_focus

    annotationScore_colors = [
        '#eaeaea', 
        apc.arcadia_all['arcadia:aster'],
        apc.arcadia_all['arcadia:aegean'],
        apc.arcadia_all['arcadia:seaweed'],
        apc.arcadia_all['arcadia:lime'],
        apc.arcadia_all['arcadia:canary']
    ]

    annotationScore_color_dict = dict(zip([str(i) for i in range(6)], annotationScore_colors))
    
    if taxon_focus == 'euk':
        taxon_color_dict = {
            'Mammalia': apc.arcadia_all['arcadia:oat'],
            'Vertebrata': apc.arcadia_all['arcadia:canary'],
            'Arthropoda': apc.arcadia_all['arcadia:seaweed'],
            'Ecdysozoa': apc.arcadia_all['arcadia:mint'],
            'Lophotrochozoa': apc.arcadia_all['arcadia:aegean'],
            'Metazoa': apc.arcadia_all['arcadia:amber'],
            'Fungi': apc.arcadia_all['arcadia:chateau'],
            'Viridiplantae': apc.arcadia_all['arcadia:lime'],
            'Sar': apc.arcadia_all['arcadia:rose'],
            'Excavata': apc.arcadia_all['arcadia:wish'],
            'Amoebazoa': apc.arcadia_all['arcadia:periwinkle'],
            'Eukaryota': apc.arcadia_all['arcadia:aster'], 
            'Bacteria': apc.arcadia_all['arcadia:slate'], 
            'Archaea': apc.arcadia_all['arcadia:dragon'],
            'Viruses': apc.arcadia_all['arcadia:orange']
        }
    elif taxon_focus == 'bac':
        taxon_color_dict = {
            'Pseudomonadota': apc.arcadia_all['arcadia:periwinkle'],
            'Nitrospirae': apc.arcadia_all['arcadia:vitalblue'],
            'Acidobacteria': apc.arcadia_all['arcadia:mars'],
            'Bacillota': apc.arcadia_all['arcadia:mint'],
            'Spirochaetes': apc.arcadia_all['arcadia:aegean'],
            'Cyanobacteria': apc.arcadia_all['arcadia:seaweed'],
            'Actinomycetota': apc.arcadia_all['arcadia:canary'],
            'Deinococcota': apc.arcadia_all['arcadia:rose'],
            'Bacteria': apc.arcadia_all['arcadia:slate'],
            'Archaea': apc.arcadia_all['arcadia:dragon'],
            'Viruses': apc.arcadia_all['arcadia:orange'],
            'Metazoa': apc.arcadia_all['arcadia:amber'],
            'Fungi': apc.arcadia_all['arcadia:chateau'],
            'Viridiplantae': apc.arcadia_all['arcadia:lime'],
            'Eukaryota': apc.arcadia_all['arcadia:wish'],
        }
    
    source_color_dict = {'blast': apc.arcadia_all['arcadia:canary'], 
                         'foldseek': apc.arcadia_all['arcadia:aegean'], 
                         'blast+foldseek': apc.arcadia_all['arcadia:amber'], 
                         'None': '#eaeaea'}
    
    plotting_rules = {
        'proteinDescription.recommendedName.fullName.value': {
            'type': 'hovertext',
            'fillna': '',
            'textlabel': 'Description'
        },
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
        'LeidenCluster': {
            'type': 'categorical',
            'fillna': 'None',
            'apply': lambda x: str(x),
            'color_order': apc.arcadia_All_ordered.values(),
            'textlabel': 'Leiden Cluster'
        },
        # 'StruCluster': {
        #     'type': 'categorical',
        #     'fillna': 'None',
        #     'color_order': apc.arcadia_All_ordered.values(),
        #     'textlabel': 'Structural Cluster'
        # },
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
        },
        'source.method': {
            'type': 'categorical',
            'fillna': 'None',
            'color_dict': source_color_dict,
            'textlabel': 'Source'
        },
    }
    
    for keyid in keyids:
        plotting_rules[f'TMscore_v_{keyid}'] = {'type': 'continuous', 'fillna': 0, 'textlabel': f'TMscore vs. {keyid}'}
        plotting_rules[f'{keyid}.hit'] = {'type': 'hovertext', 'apply': lambda x: True if x == 1 else False, 'fillna': 'None', 'textlabel': f'Blast/Foldseek Hit to {keyid}'}
    
    coordinates_file = apply_coordinates(dimensions_file, features_file, save = True, prep_step = True, 
                                         dimtype = dimensions_type)
    plot_interactive(coordinates_file, plotting_rules, output_file = output_file, keyids = keyids)
    
if __name__ == '__main__':
    main()
