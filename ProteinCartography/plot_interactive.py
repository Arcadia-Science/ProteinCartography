#!/usr/bin/env python
import argparse
import plotly
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import arcadia_pycolor as apc
import numpy as np
import matplotlib.colors as mc
import colorsys

# only import these functions when using import *
__all__ = ["adjust_lightness", "apply_coordinates", "assign_taxon", "extend_colors", "plot_interactive"]

# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dimensions", required = True, help = "path to dimensionality reduction file")
    parser.add_argument("-f", "--features", required = True, help = "path to aggregated features file")
    parser.add_argument("-o", "--output", required = True, help = "path of destination HTML file")
    parser.add_argument("-t", "--dimensions-type", default = '', help = "dimensions data type (pca, tsne, umap, etc.)")
    parser.add_argument("-k", "--keyids", nargs = "+", default = [], help = "keyids for plotting. Usually input protids.")
    parser.add_argument("-x", "--taxon_focus", default = "euk", help = "Coloring scheme/ taxonomic groups for broad taxon plot.\n'euk'(aryote) is default; 'bac'(teria) is another option.")
    parser.add_argument("-X", "--plot-width", default = '700', help = "width of resulting plot.")
    parser.add_argument("-Y", "--plot-height", default = '700', help = "width of resulting plot.")
    args = parser.parse_args()
    
    return args

def adjust_lightness(color: str, amount = 0.5) -> str:
    '''
    Takes a HEX code or matplotlib color name and adjusts the lightness.
    Values < 1 result in a darker color, whereas values > 1 result in a lighter color.
    
    Args:
        color (str): hex value (e.g. "#FEACAF") or a valid matplotlib color name (e.g. "tab:blue").
        amount (float): values < 1 produce a darker color and values > 1 produce a lighter color.
    Returns:
        resulting color as HEX string.
    '''
    try:
        color_string = mc.cnames[color]
    except:
        color_string = color
    # convert rgb to hls
    color_string = colorsys.rgb_to_hls(*mc.to_rgb(color_string))
    # adjust the lightness in hls space and convert back to rgb
    color_string2 = colorsys.hls_to_rgb(color_string[0], max(0, min(1, amount * color_string[1])), color_string[2])
    # return the new rgb as a hex value
    return mc.to_hex(color_string2)

def apply_coordinates(dimensions_file: str, features_file: str, saveprefix = None, dimtype = None, save = False, prep_step = False):
    '''
    Adds coordinate information for data from the dimensions file to the aggregated features file.
    
    Args:
        dimensions_file (str): path to dimensions file.
        features_file (str): path to aggregated features file.
        saveprefix (str): prefix for output file.
        dimtype (str): type of dimensions file (e.g. pca, tsne, etc.). If not provided, is inferred from dimensions file.
        save (bool): whether or not to save the file. Defaults to False.
        prep_step (bool): if True, returns path of the final output file. Defaults to False.
    '''
    
    # Load the dimensions data
    reduced_dim_df = pd.read_csv(dimensions_file, sep = '\t')
    
    # If dimtype not provided, infer based on second column of dataframe
    if dimtype is None:
        outdimtype = ''.join([i for i in reduced_dim_df.columns[1] if not i.isdigit()])
    else:
        outdimtype = dimtype
    
    # Check to make sure there's a tsv there
    if '.tsv' not in features_file:
        raise Exception(f'{features_file} does not end in ".tsv" as expected.')
    
    # Load the features data
    agg_features_df = pd.read_csv(features_file, sep = '\t')
    
    # Merge the two dataframes
    plot_features_df = reduced_dim_df.merge(agg_features_df, on = 'protid')
    
    # Generate filepath based on dimensions type and prefix
    if saveprefix is not None:
        savefile = '_'.join([saveprefix, outdimtype + '.tsv'])
    else:        
        savefile = features_file.replace('.tsv', '_' + outdimtype + '.tsv')
    
    # Save if needed
    if save:
        plot_features_df.to_csv(savefile, sep = '\t', index = None)
    
    # If prep_step, return the savefile path
    # Otherwise, return the results
    if prep_step:
        return savefile
    else:
        return plot_features_df

def assign_taxon(taxon_list: list, rank_list: list, hierarchical = False, sep = ','):
    '''
    Takes a list of taxa and assigns it a category based on a ranked-order list.
    
    The ranked-order list should be ordered from most-specific to least-specific taxonomic groupings.
    The function iterates through the list of ranks and checks if each is within the taxon list.
    If hierarchical = True, returns only the first element. Otherwise, concatenates the hits in order.
    
    Args:
        taxon_list (list): list of taxa for the object in question.
        rank_list (list): ranked list to search through.
        hierarchical (bool): whether to return the first element found in the rank_list or a concatentation of all hits.
        sep (str): separator character when returning a concatenation of hits.
    Returns:
        a string either containing the first-rank hit or the string concatenation of all hits.
    '''
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
    
def extend_colors(color_keys: list, color_order: list, steps = [0.7, 0.5, 0.3]) -> list:
    '''
    Checks a list of keys and colors and extends the colors list as needed.
    
    Args:
        color_keys (list): list of color keys
        color_order (list): list of HEX colors
        steps (list): float list of potential lightness adjustments to make
    Returns:
        extended list of colors, not including original colors
    '''
    num_cycles = int(np.round(len(color_keys) / len(color_order)))

    # collector for additional colors
    more_colors = []

    # if there aren't enough cycles, die
    if num_cycles > len(steps):
        raise Exception(f'Can create up to {len(steps) * len(color_order)} colors to use.\nNeeded {len(color_keys)} colors.')

    # create additional colors and add to collector
    for n in range(num_cycles - 1):
        color_order_duplicated = [adjust_lightness(color) for color in color_order]
        more_colors.extend(color_order_duplicated)
        
    return more_colors

def plot_interactive(coordinates_file: str, plotting_rules: dict,
                    marker_size = 4, marker_opacity = 0.8, output_file = '', keyids = [], show = False,
                    plot_width = 700, plot_height = 700):
    '''
    Plots all proteins on a 2D interactive Plotly plot using a set of rules.
    
    The plotting_rules dictionary should have the following format.
    Each column is an entry in the dictionary containing a dictionary of rules.
    ```
    {
        'column1.name': {
            'type': 'categorical',
            'parameter1': value,
            'parameter2': value,
            ...
        }
        'column2.name': {
            'type': 'hovertext',
            ...
        }
    }
    ```
    The possible rules for each column are as follows:
    ### For any plot type ###
    - **'type' (required):**
        - `'categorical'`, `'continuous'`, `'taxonomic'`, or `'hovertext'`
        - Determines the plotting style of the data.
        - If 'categorical', expects the data to be in string format.
        - If 'continuous', expects the data to be in numerical (int or float) format.
        - If 'taxonomic', expects an ordered list of taxa to be used for grouping.
        - If 'hovertext', won't plot the data as a dropdown menu option, but will include the data in the hover-over text box.
    - **'fillna':**
        - A value to fill rows that are `np.nan`.
        - For categorical plots, usually empty string `''`.
        - For continuous plots, usually `0`.
    - **'apply':**
        - A function that will be applied to every element of the feature before plotting.
        - This can be used to convert a numerical value into a categorical value with a lambda function.
        - For example, `lambda x: str(x)`
    - **'skip_hover':**
        - Boolean, whether to include this data in the hovertext.
    - **'textlabel':**
        - A string that replaces the column name on buttons and in hover text.
        - Useful if the column name is too ugly.
        
    ### For 'categorical' and 'taxonomic' plots ###
    - **'color_order':**
        - A list of HEX code colors used for coloring categorical variables.
        - If there aren't enough colors for unique values of the data, will generate up to 3x more colors of varying darkness.
    - **'color_dict':**
        - A dictionary of key: value pairs for coloring categorical variables.
        - The key is the name of the category and the value is a HEX code color.
    
    ### For 'taxonomic' plots ###
    - **'taxon_order':**
        - Exclusively used for 'taxonomic' style plots. A list of ranked-order taxa for categorization.
        - Should start with more-specific taxa and expand to less-specific taxa.
    
    Args:
        coordinates_file (str): path to coordinates file (an aggregated features file with columns for the coordinates in PC/TSNE/UMAP space).
        plotting_rules (dict): a JSON-like dictionary containing plotting rules for relevant columns of the data.
        marker_size (int): size of markers.
        marker_opacity (float): opacity of markers (value should between 0 and 1).
        output_file (str): path of destination file.
        keyids (list): list of key protids to assign a star-shaped marker. Usually input proteins.
        show (bool): whether or not to show the plot.
    Returns:
        if show = False, returns the plotly.graphobjects object of the plot.
    '''
    # make a copy of the starting data
    df = pd.read_csv(coordinates_file, sep = '\t')
    
    # generate a hover template text string
    # custom data in plotly is indexed by order
    hovertemplate_generator = ["<b>%{customdata[0]}</b></br>––––––––––––"]
    # pass in the custom data columns
    # this sets protid as customdata[0]
    # and the rest of the plotting rules in order
    
    valid_plotting_rules = [col for col in plotting_rules.keys() if col in df.columns]
    custom_data = ['protid'] + list(valid_plotting_rules)
    
    # iterate through plotting rules, applying relevant rules
    for col in plotting_rules:
        if col not in df.columns:
            continue
        
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
        if col not in df.columns:
            continue
        
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
                # if not, extend available colors
                if len(color_order) < len(color_keys):
                    more_colors = extend_colors(color_keys, color_order, steps = [0.7, 0.5, 0.3])
                    
                    # extend color order with new colors
                    color_order.extend(more_colors)
                
                # make color dict
                colors_dict = dict(zip(color_keys, color_order))
            
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
        
        # Plotting rules for taxonomic plots
        elif plotting_rules[col]['type'] == 'taxonomic':
            
            # Check to make sure there's a taxon order
            if 'taxon_order' in plotting_rules[col].keys():
                taxon_order = plotting_rules[col]['taxon_order']
            else:
                raise Exception('Please provide a "taxon_order" list in the plotting rules.')
            
            # Set color keys to taxon order
            color_keys = taxon_order
            
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
                # if not, extend available colors
                if len(color_order) < len(color_keys):
                    more_colors = extend_colors(color_keys, color_order, steps = [0.7, 0.5, 0.3])
                    
                    # extend color order with new colors
                    color_order.extend(more_colors)
                
                # make colors dict
                colors_dict = dict(zip(color_keys, color_order))
                
                # add entry for 'Other'
                colors_dict['Other'] = '#eaeaea'
            
            # generate a taxonomic category column using assign_taxon lambda function
            col_taxonomic = col + '_taxonomic'
            df[col_taxonomic] = df[col].apply(lambda x: assign_taxon(x, taxon_order, hierarchical = True))
            
            # plot using the above parameters
            plots[col] = px.scatter(df, dim1, dim2, color = col_taxonomic, hover_name = 'protid',
                                   color_discrete_map = colors_dict,
                                   custom_data = custom_data)
            plots[col].update_traces(marker = dict(size = marker_size, opacity = marker_opacity),
                                    hovertemplate = hovertemplate)
    
    # collector for the number of scatter items
    scatter_counter = {}
    # collector for the order of figures added to the plot
    fig_order = []
    
    # create new empty figure to move every original figure onto
    # this is a workaround to allow colored legends to be preserved and automatically switched using the dropdown
    # the alternative approach forces you to show all legend items at the same time for all colorations, which is messy
    fig = go.Figure()
    
    # iterate through the plots to get their objects
    for j, (col, plot) in enumerate(plots.items()):
        
        # if it's the first plot, make sure it's visible
        if j == 0:
            vis = True
        # otherwise, hide it
        else:
            vis = False
        
        # appends the figure's name to the list of figure orders
        fig_order.append(col)
        
        # counts the number of different scatter traces within the original plot
        # this is needed to determine the number of true or false values for the dropdown toggle visibility field
        scatter_counter[col] = len(plot.data)
        
        # if there are more than 15 different categories, hide the legend
        if scatter_counter[col] > 15:
            showlegend = False
        else:
            showlegend = True
        
        # finally, add the original plot to the new plot.
        for scatter in plot.data:
            
            if type(scatter) == plotly.graph_objs._scattergl.Scattergl:
                fig.add_trace(go.Scattergl(scatter, visible = vis, showlegend = showlegend))
            elif type(scatter) == plotly.graph_objs._scatter.Scatter:
                fig.add_trace(go.Scatter(scatter, visible = vis, showlegend = showlegend))
    
    # if there are any keyids provided, generate an additional plot
    if keyids != []:
        
        # get the positions of the key ids
        keypoints = df[df['protid'].isin(keyids)]
        
        # add those key points as a new scatter with a star-diamond marker
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
    
    # define a helper function to determine the visibility toggle list
    # this sets, for every single trace (each category of each plot gets a trace), its visibility
    def visibility_list(col):
        entries = []
        
        for fig in fig_order:
            if fig == col:
                entries.extend([True] * scatter_counter[fig])
            else:
                entries.extend([False] * scatter_counter[fig])
        
        return entries
    
    # collector for buttons
    buttons = []
    
    # for each plot, create a button
    for k, (col, plot) in enumerate(plots.items()):
        
        # get the textlabel for that plot if it's provided
        # otherwise, default to column name
        if 'textlabel' in plotting_rules[col].keys():
            button_label = plotting_rules[col]['textlabel']
        else:
            button_label = col
        
        # create a button that toggles visibility based on which data is currently selected
        button_item = dict(
                        args = ["visible", visibility_list(col) + [True]],
                        label = button_label,
                        method = "restyle"
                    )
        # add this button to the selector
        buttons.append(button_item)
        
    # collector for dropdown menu
    dropdown_menu = dict(
                buttons = list(buttons),
                showactive = True,
                x = 0.1,
                xanchor = "left",
                y = 1.1,
                yanchor = "top"
    )
    
    # create a separate button to toggle keyid trace
    # the visibility will match the visibility of the first plot upon generation
    # I couldn't find an obvious way to make ONLY this plot's visibility toggle while preserving the rest...
    # current behavior is when this button is clicked, the visibility of the keyid stars is toggled
    # but all other data except the first coloration style are hidden
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
        # add this to the menu buttons
        updatemenu = [dropdown_menu, keyids_button]
    else:
        # if no keyids, don't add this button
        updatemenu = [dropdown_menu]
    
    # add the menu buttons
    fig.update_layout(
        updatemenus = updatemenu
    )
    
    # set the width, height, and other parameters of the plot
    fig.update_layout(
        width = plot_width, height = plot_height,
        annotations=[
            dict(text="color", x=0.01, xref="paper", y=1.08, yref="paper",
                 align="left", showarrow=False) # This shows the word "color" next to the dropdown
        ],
        plot_bgcolor= '#fafafa'
    )
    
    # hide tick labels
    fig.update_layout(xaxis = dict(showticklabels = False), yaxis = dict(showticklabels = False))
    # make x and y axes have a fixed scale
    fig.update_yaxes(scaleanchor = "x", scaleratio = 1)
    
    # move the legend for categorical and color bar features to the bottom
    fig.update_layout(legend=dict(orientation = "h", yanchor="top", y = 0, xanchor="left", x = 0, font = dict(size = 10)))
    fig.update_layout(coloraxis_colorbar=dict(yanchor="top", y=0, x=0.5,
                                          ticks="outside", orientation = 'h'))
    
    # save if filename is provided
    if output_file != '':
        fig.write_html(output_file)
    
    # show if desired
    if show:
        fig.show()
    
    # return figure object with everything in it
    return fig

# run this if called from the interpreter
def main():
    args = parse_args()
    dimensions_file = args.dimensions
    features_file = args.features
    output_file = args.output
    dimensions_type = args.dimensions_type
    keyids = args.keyids
    taxon_focus = args.taxon_focus
    plot_width = int(args.plot_width)
    plot_height = int(args.plot_height)
    
    # color order for annotation score
    annotationScore_colors = [
        '#eaeaea', 
        apc.arcadia_all['arcadia:aster'],
        apc.arcadia_all['arcadia:aegean'],
        apc.arcadia_all['arcadia:seaweed'],
        apc.arcadia_all['arcadia:lime'],
        apc.arcadia_all['arcadia:canary']
    ]
    
    # color dictionary for annotation score (values from 0 to 5)
    annotationScore_color_dict = dict(zip([str(i) for i in range(6)], annotationScore_colors))
    
    # if the taxonomic focus is eukaryote, use these groupings and colors
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
            'Viruses': apc.arcadia_all['arcadia:denim']
        }
        
    # else if the taxonomic focus is bacteria, use these groupings and colors
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
            'Viruses': apc.arcadia_all['arcadia:denim'],
            'Metazoa': apc.arcadia_all['arcadia:amber'],
            'Fungi': apc.arcadia_all['arcadia:chateau'],
            'Viridiplantae': apc.arcadia_all['arcadia:lime'],
            'Eukaryota': apc.arcadia_all['arcadia:wish'],
        }
    
    # use this color dictionary for the sources (blast vs foldseek vs blast+foldseek)
    source_color_dict = {'blast': apc.arcadia_all['arcadia:canary'], 
                         'foldseek': apc.arcadia_all['arcadia:aegean'], 
                         'blast+foldseek': apc.arcadia_all['arcadia:amber'], 
                         'None': '#eaeaea'}
    
    # use this dictionary as default plotting rules
    # the order of elements in this dictionary determines their order 
    # both in the plotting dropdown list and the hovertext
    plotting_rules = {
        'proteinDescription.recommendedName.fullName.value': {
            'type': 'hovertext', # this data only shows up in hovertext
            'fillna': '',
            'textlabel': 'Description'
        },
        'organism.scientificName': {
            'type': 'hovertext', # this data only shows up in hovertext
            'fillna': '',
            'textlabel': 'Species'
        },
        'organism.commonName': {
            'type': 'hovertext', # this data only shows up in hovertext
            'fillna': '',
            'textlabel': 'Common name'
        },
        'LeidenCluster': {
            'type': 'categorical',
            'fillna': 'None',
            'apply': lambda x: str(x), # Leiden cluster is often read as int; this forces it to be string
            'color_order': apc.arcadia_All_ordered.values(),
            'textlabel': 'Leiden Cluster'
        },
        # 'StruCluster': { # This is hidden for now as people thought it wasn't very useful
        #     'type': 'categorical',
        #     'fillna': 'None',
        #     'color_order': apc.arcadia_All_ordered.values(),
        #     'textlabel': 'Structural Cluster'
        # },
        'annotationScore': {
            'type': 'categorical',
            'fillna': 0,
            'apply': lambda x: str(int(x)), # Annotation score is parsed as float but we want it to be string
            'color_dict': annotationScore_color_dict,
            'textlabel': 'Annotation Score'
        },
        'organism.lineage': {
            'type': 'taxonomic',
            'fillna': '[]',
            'apply': lambda x: eval(x), # This converts the taxonomic groupings from a string-ified list to a real list
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
    
    # Add additional plotting rules for each keyid (usually the input protein)
    for keyid in keyids:
        # Add a plot for tmscore to each input protein
        plotting_rules[f'TMscore_v_{keyid}'] = {
            'type': 'continuous', 
            'fillna': 0, 
            'textlabel': f'TMscore vs. {keyid}'
        }
        # Add hovertext for whether or not a given protein was a hit via blast or foldseek to the input protein
        plotting_rules[f'{keyid}.hit'] = {
            'type': 'hovertext', 
            'apply': lambda x: True if x == 1 else False, 
            'fillna': 'None', 
            'textlabel': f'Blast/Foldseek Hit to {keyid}'
        }
    
    # generate coordinates file for the plot
    coordinates_file = apply_coordinates(dimensions_file, features_file, save = True, prep_step = True, 
                                         dimtype = dimensions_type)
    # make the plot
    plot_interactive(coordinates_file, plotting_rules, 
                     output_file = output_file, keyids = keyids,
                     plot_width = plot_width, plot_height = plot_height)

# check if called from interpreter
if __name__ == '__main__':
    main()
