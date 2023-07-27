#!/usr/bin/env python
import argparse
from dash import Dash, html, dash_table, dcc, Input, Output, State, callback, ctx
import dash_bio as dashbio
import dash_daq as daq
from dash_bio.utils import PdbParser, create_mol3d_style
import pandas as pd
import numpy as np
import arcadia_pycolor as apc
import matplotlib.colors as mcolors
import colorsys
import json
from pathlib import Path
import os

from plot_interactive import generate_plotting_rules, plot_interactive, assign_taxon
from cluster_similarity import plot_group_similarity
from pdb_tools import assign_residue_colors, extract_residue_confidence, RESIDUE_CONFIDENCE_COLORS

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-path", required = True, help = "path to ProteinCartography output folder")
    parser.add_argument("-n", "--analysis-name", default = '', help = "name of the analysis")
    parser.add_argument("-t", "--dimensions-type", default = '', help = "dimensions data type (pca, tsne, umap, etc.)")
    parser.add_argument("-k", "--keyids", nargs = "+", default = [], help = "keyids for plotting. Usually input protids.")
    parser.add_argument("-x", "--taxon_focus", default = "euk", help = "Coloring scheme/ taxonomic groups for broad taxon plot.\n'euk'(aryote) is default; 'bac'(teria) is another option.")
    args = parser.parse_args()
    
    return args

def make_data_dict(input_path: str, 
                   analysis_name = '', taxon_focus = 'euk', keyids = [],
                   plot_width = 500, plot_height = 600):
    
    input_Path = Path(input_path) 

    results_Path = input_Path / 'clusteringresults'

    if analysis_name == '':
        analysis_name = [i for i in os.listdir(results_Path) if '_aggregated_features.tsv' in i][0].replace('_aggregated_features.tsv', '')

    tsne_file = results_Path / f'{analysis_name}_aggregated_features_pca_tsne.tsv'
    umap_file = results_Path / f'{analysis_name}_aggregated_features_pca_umap.tsv'
    similarity_file = results_Path / f'{analysis_name}_leiden_similarity.tsv'

    plotting_rules = generate_plotting_rules(taxon_focus, keyids)

    structures_Path = input_Path / 'foldseekclustering'

    if os.path.exists(structures_Path):
        structures_list = [file for file in os.listdir(structures_Path)]
    else:
        structures_list = []

    scatter_tsne = plot_interactive(tsne_file, plotting_rules, keyids = keyids,
                    plot_width = plot_width, plot_height = plot_height,
                    show = False, marker_size = 4, marker_opacity = 0.8,
                    plot_bgcolor='white', hide_hover = True)
    
    scatter_umap = plot_interactive(umap_file, plotting_rules, keyids = keyids,
                plot_width = plot_width, plot_height = plot_height,
                show = False, marker_size = 4, marker_opacity = 0.8,
                plot_bgcolor='white', hide_hover = True)
    
    scatter_data = pd.read_csv(tsne_file, sep = '\t').drop_duplicates()

    heatmap_figure = plot_group_similarity(
        similarity_file, show = False, plot_width = 450, plot_height = 450
        )
    
    data_dict = {
        'analysis_name': analysis_name,
        'scatterplot_tsne': scatter_tsne,
        'scatterplot_umap': scatter_umap,
        'scatter_data': scatter_data,
        'plotting_rules': plotting_rules,
        'heatmap': heatmap_figure,
        'structures_path': structures_Path,
        'structures_list': structures_list,
        'keyids': keyids
    }

    return data_dict

def css_linear_gradient_cutoff(tuple_list: list, max_bound_percentage: int, max_bound_color = '#FFFFFF'):
    percent_list = [(color, np.round(value * 100, 0)) for (value, color) in tuple_list]
    cutoff_list = {color: [value] for (color, value) in percent_list if value < max_bound_percentage}
    
    last_entry = list(cutoff_list.keys())[-1]
    cutoff_list[last_entry].append(max_bound_percentage)
    cutoff_list[max_bound_color] = [max_bound_percentage]
    
    css_string = ', '.join([' '.join([color, ' '.join([f'{val}%' for val in values])]) for color, values in cutoff_list.items()])
    
    return css_string

def data_bars(df, column, grad_tuple_list):
    n_bins = 100
    bounds = [i * (1.0 / n_bins) for i in range(n_bins + 1)]
    ranges = [
        ((df[column].max() - df[column].min()) * i) + df[column].min()
        for i in bounds
    ]
    styles = []
    for i in range(1, len(bounds)):
        min_bound = ranges[i - 1]
        max_bound = ranges[i]
        max_bound_percentage = bounds[i] * 100
        styles.append({
            'if': {
                'filter_query': (
                    '{{{column}}} >= {min_bound}' +
                    (' && {{{column}}} < {max_bound}' if (i < len(bounds) - 1) else '')
                ).format(column=column, min_bound=min_bound, max_bound=max_bound),
                'column_id': column
            },
            'background': (
                f'linear-gradient(90deg, {css_linear_gradient_cutoff(grad_tuple_list, max_bound_percentage)})'
            ),
            'paddingBottom': 2,
            'paddingTop': 2
        })

    return styles

def search_df(df: pd.DataFrame, query: str, low_memory = False):
    '''
    Takes an input TSV filepath and searches for any rows that contain the query string.
    
    Args:
        tsv (str): path to input TSV file to search through.
        query (str): substring to search for in all columns and rows.
    Returns:
        pandas.DataFrame containing rows that matched the query.
    '''
    searched_df = df[df.apply(lambda r: r.str.contains(query, case=False).any(), axis=1)]
    
    return searched_df

def make_dash_app(data_dict: dict):
    # take input data
    analysis_name = data_dict['analysis_name']
    scatter_tsne = data_dict['scatterplot_tsne']
    scatter_umap = data_dict['scatterplot_umap']
    scatter_data = data_dict['scatter_data']
    plotting_rules = data_dict['plotting_rules']
    heatmap_figure = data_dict['heatmap']
    structures_path = data_dict['structures_path']
    structures_list = data_dict['structures_list']
    keyids = data_dict['keyids']

    # Initialize the app
    app = Dash(__name__)

    scatter_data['LeidenCluster_circle'] = "● " + scatter_data['LeidenCluster']
    scatter_data['Annotation_circle'] = "● " + scatter_data['Annotation'].astype(int).astype(str)

    taxon_order = plotting_rules['Lineage']['taxon_order']
    taxon_color_dict = dict(zip(taxon_order, plotting_rules['Lineage']['color_order']))
    taxon_color_dict['Other'] = apc.All['arcadia:brightgrey']
    scatter_data['Taxon'] = scatter_data['Lineage'].apply(lambda x: assign_taxon(x, taxon_order, hierarchical = True))
    scatter_data['Taxon_circle'] = "● " + scatter_data['Taxon']

    source_color_dict = plotting_rules['source.method']['color_dict']
    
    uniprot_set = set(scatter_data['Entry'].values)

    def link_generator(row):
        protid = row['protid']
        if row['Entry'] is not np.nan:
            return f'[{protid}](https://www.uniprot.org/uniprotkb/{protid})'
        else:
            return 'Not in Uniprot'
        
    scatter_data['protid_link'] = scatter_data.apply(lambda row: link_generator(row), axis = 1)

    column_name_mapper = {
        'Protein names': 'Protein name',
        'Gene Names (primary)': 'Gene name',
        'Organism': 'Organism',
        'LeidenCluster_circle': 'Cluster',
        'Annotation_circle': 'Score',
        'Taxon_circle': 'Taxon',
        'Length': 'Length'
    }

    default_columns = list(column_name_mapper.keys())

    protid_columns = [{'id': 'protid', 'name': 'protid',}, {'id': 'protid_link', 'name': 'Uniprot', 'type': 'text', 'presentation': 'markdown'}]

    data_table_columns = protid_columns + [{'id': c, 'name': column_name_mapper[c]} for c in default_columns]

    LC_conditional_color_dict = dict(zip(sorted(
        scatter_data['LeidenCluster'].unique()), 
        plotting_rules['LeidenCluster']['color_order']
        ))
    
    LC_conditional_format = [{
            'if': {
                'filter_query': f'{{LeidenCluster}} = {key}',  # matching rows of a hidden column with the id, `id`
                'column_id': 'LeidenCluster_circle'
            },
            'color': value,
            'font-weight': 'bold'
        } for key, value in LC_conditional_color_dict.items()]
    
    Annotation_conditional_color_dict = plotting_rules['Annotation']['color_dict'].items()
    
    Annotation_conditional_format = [{
        'if': {
            'filter_query': f'{{Annotation}} = {key}',  # matching rows of a hidden column with the id, `id`
            'column_id': 'Annotation_circle'
        },
        'color': value,
        'font-weight': 'bold'
    } for key, value in Annotation_conditional_color_dict]

    Taxon_conditional_format = [{
        'if': {
            'filter_query': f'{{Taxon}} = {key}',  # matching rows of a hidden column with the id, `id`
            'column_id': 'Taxon_circle'
        },
        'color': value,
        'font-weight': 'bold'
    } for key, value in taxon_color_dict.items()]
    
    Length_color_scale = apc.Gradients['arcadia:cividis_r'].mpl_LinearSegmentedColormap.resampled(10)
    Length_color_list = [(i / Length_color_scale.N, mcolors.rgb2hex(Length_color_scale(i))) for i in range(Length_color_scale.N)]
    
    Length_conditional_format = data_bars(scatter_data, 'Length', Length_color_list)

    cell_conditional_formats = LC_conditional_format + Annotation_conditional_format + Taxon_conditional_format + Length_conditional_format 

    tooltip_maker = lambda x: [
                        { column: {'value': str(value), 'type': 'markdown'}
                        for column, value in row.items()
                        } for row in x.to_dict('records')
                    ]
    tooltip_data = tooltip_maker(scatter_data)
    scatter_records = scatter_data.to_dict('records')

    cell_style =  {
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                        'maxWidth': 0,
                        'textAlign': 'left',
                        'font-family': 'Arial',
                        'padding-right': '5px',
                        'padding-left': '5px',
                        'padding-top' : '2px',
                        'padding-bottom' : '2px',
                    }

    scatter_config = {'displayModeBar': True,
                      'toImageButtonOptions': {
                          'filename': f'{analysis_name}_scatter',
                          'format': 'svg'
                          },
                        'modeBarButtonsToRemove': ['zoomIn', 'zoomOut'],
                        'displaylogo' : False,
                    }
    
    heatmap_config = {'displayModeBar': False}
    
    tabs_styles = {
        'height': '40px',
        'marginTop': '-14px',
        'marginLeft': '-12px',
        'marginRight': '-12px'
    }
    tab_style = {
        'fontFamily': 'Arial, sans-serif',
        'color': '#B6C8D4',
        'backgroundColor': '#FEFEFE',
        'fontSize': '20px',
        'fontWeight': 'lighter',
        'letterSpacing': '2px',
        'textAlign': 'center',
        'verticalAlign': 'middle',
        'paddingTop': '8px',
        'paddingLeft': '0px',
        'paddingRight': '0px',
        'paddingBottom': '0px',
        'border': '0px',
        'boxShadow': '1px 1px 2px 1px #EAEAEA' 
    }
    tab_selected_style = {
        'fontFamily': 'Arial, sans-serif',
        'color': '#000000',
        'fontSize': '20px',
        'fontWeight': 'lighter',
        'letterSpacing': '2px',
        'border': '0px',
        'borderTop': '2px solid #5088C5',
        'backgroundColor': '#FCFCFC',
        'textAlign': 'center',
        'verticalAlign': 'middle',
        'paddingTop': '8px',
        'paddingLeft': '0px',
        'paddingRight': '0px',
        'paddingBottom': '0px'
    }

    first_input_protein = keyids[0] + '.pdb'
    residue_confidence_key = dict(zip(['Very high (>90)', 'Confident (>70)', 'Low (>50)', 'Very low (<50)'], RESIDUE_CONFIDENCE_COLORS.values()))

    if structures_path is not None:
        structure_tab_content = [
                    html.H3(children = 'Structure viewer'),
                    html.Div(
                        id = 'structure-tab-menu',
                        children = [
                            dcc.Dropdown(
                                id = 'structure-protein-1-input',
                                value = first_input_protein,
                                options = structures_list,
                                maxHeight = 200
                            ),
                            #(
                            #     id = 'structure-protein-2',
                            #     type = 'search', 
                            #     children = 'Protein 2',
                            #     list = structures_list
                            # ),
                            html.Button(
                                id = 'structure-display-button',
                                children = 'Display'
                            ),
                    ]),
                    dashbio.Molecule3dViewer(
                        id = 'dashbio-mol3dviewer',
                        modelData = {'atoms': [], 'bonds': []},
                        zoom = {'factor': 1.3},
                        height = 400
                    ),
                    html.Div(
                        id = 'structure-quality-key',
                        children = [html.A('pLDDT', id = 'plddt-link', href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3799472/", target="_blank")] + [
                            html.Div([
                                html.Div(
                                    style = {'backgroundColor': color, 'width': 20, 'height': 20}), 
                                html.P(quality)
                        ]) for quality, color in residue_confidence_key.items()
                        ])
                ]
    else:
        structure_tab_content = [
                    html.H3(children = 'Structure viewer'),
                    html.P(children = 'No data provided.')
                ]

    # App layout
    app.layout = html.Div([
        html.Header(
            id = 'title-bar',
            children = [
                html.Div(
                    id = 'title-bar-left',
                    children = [
                    html.H1(
                        id = 'logo',
                        children = 'ProteinCartography',
                    ),
                    html.P(
                        id = 'title-tag',
                        children = 'Analysis',
                    ),
                    html.Code(
                        id = 'title-code',
                        children = analysis_name
                    )
                ]),
                html.Div(
                    id = 'title-bar-right',
                    children = [
                        html.A(
                            id = 'github-a',
                            href="https://www.github.com/Arcadia-Science/gene-family-cartography",
                            target = "_blank",
                            children = [
                                html.Code(id = 'version', children = 'v0.3.0'), 
                                html.Img(id = 'github-logo', src = 'assets/github-mark.svg')
                            ]
                        ),
                        html.A(
                            id = 'arcadia-a',
                            href="https://www.arcadiascience.com/",
                            target = "_blank",
                            children = [
                                html.Img(id = 'arcadia-logo', src = 'assets/arcadia-logo-white.svg')
                            ]
                        ),
                ])
        ]),
        html.Div(
            id = 'top-card-row',
            children = [
                html.Div(
                    id = 'info-card',
                    children = [
                        html.Div(
                            id = 'info-card-titlebar',
                            children = [
                                html.H2(
                                    id = 'info-card-title',
                                    children = 'Info'
                                ),
                                html.Button(
                                    id = 'tutorial-button',
                                    children = 'Tutorial',
                                    n_clicks = 0
                                )   
                        ]),
                        dcc.Markdown(
                            id = 'info-text',
                            children = '',
                            style = {'white-space': 'pre-wrap'},
                            dangerously_allow_html = True
                        )
                    ]
                ),
                html.Div(
                    id = 'scatter-card',
                    children = [
                        html.Div(
                            id = 'scatter-card-titlebar',
                            children = [
                                html.H2(
                                    id = 'scatter-card-title',
                                    children = 'Scatter'
                                ),
                                html.Div(
                                    id = 'space-toggle',
                                    children = [
                                        html.A('t-SNE', href="https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html", target="_blank", style = {'marginRight': '5px'}),
                                        daq.ToggleSwitch(
                                            id = 'space-toggle-switch',
                                            value = False,
                                            size = 30,
                                        ),
                                        html.A('UMAP', href="https://umap-learn.readthedocs.io/en/latest/", target="_blank", style = {'marginLeft': '5px'})
                                    ]
                                )   
                            ]
                        ),
                        dcc.Graph(
                            id = 'scatter-plot',
                            figure = scatter_tsne,
                            config = scatter_config
                        ),
                ]),
                html.Div(
                    id = 'analysis-card',
                    children = [
                        dcc.Tabs(
                            id = 'analysis-tabs',
                            value = 'heatmap-tab',
                            style = tabs_styles,
                            children = [
                                dcc.Tab(
                                    label = 'Heatmap', 
                                    value = 'heatmap-tab', 
                                    style = tab_style, 
                                    selected_style = tab_selected_style,
                                    children = [
                                        html.Div(
                                            id = 'heatmap-tab-content',
                                            children = [
                                                html.H3(children = 'Cluster similarity heatmap'),
                                                dcc.Graph(figure = heatmap_figure, config = heatmap_config)
                                            ])
                                ]),
                                dcc.Tab(
                                    label = 'Structure', 
                                    value = 'structure-tab', 
                                    style = tab_style, 
                                    selected_style = tab_selected_style,
                                    children = [
                                        html.Div(
                                            id = 'structure-tab-content',
                                            children = structure_tab_content
                                        )
                                ]),
                                dcc.Tab(
                                    label = 'Semantics', 
                                    value = 'semantics-tab', 
                                    style = tab_style, 
                                    selected_style = tab_selected_style,
                                    children = [
                                        html.Div(
                                            id = 'semantics-tab-content',
                                            children = [
                                                html.H3(children = 'Semantic analysis'),
                                                html.P(children = 'In development :)')
                                            ])
                                ]),
                                dcc.Tab(
                                    label = 'Phylogeny', 
                                    value = 'phylogeny-tab', 
                                    style = tab_style, 
                                    selected_style = tab_selected_style,
                                    children = [
                                        html.Div(
                                            id = 'phylogeny-tab-content',
                                            children = [
                                                html.H3(children = 'Phylogeny'),
                                                html.P(children = 'In development :)')
                                            ])
                            ])
                        ])
                ]),
        ]),
        html.Div(
            id = 'database-card', 
            children = [
                html.Div(
                    id = 'database-card-titlebar',
                    children = [
                        html.H2(
                            id = 'database-card-title',
                            children = 'Database'
                        ),
                        html.Div(
                            id = 'search-div',
                            children = [
                                dcc.Input(
                                    id = 'search-bar', type = 'text', value = '', n_submit = 0,
                                ),
                                html.Button(
                                    id = 'search-button', n_clicks = 0, children = 'Search'
                                ),
                                html.Button(
                                    id = 'reset-button', n_clicks = 0, children = 'Reset'
                                ),
                                html.Div(
                                    id = 'download-div',
                                    children = [
                                        html.Button(
                                            id = 'download-button', n_clicks = 0, children = 'Download'
                                        ),
                                        dcc.Download(
                                            id = 'download-tsv-component'
                                    )],
                                )

                        ])
                ]),
                dash_table.DataTable(
                    id = 'data-table',
                    data = scatter_records,
                    columns = data_table_columns,
                    page_size = 10,
                    style_header = {
                        'font-size': 14,
                        'font-weight': 'bold'
                    },
                    style_data = {
                        'whitespace': 'normal',
                        'font-size': 12
                    },
                    style_cell = cell_style,
                    style_data_conditional = cell_conditional_formats,
                    style_cell_conditional = [
                        {
                            'if': {'column_id': 'protid'},
                            'textAlign': 'right',
                            'width' : '8%'
                        },
                        {
                            'if': {'column_id': 'protid_link'},
                            'textAlign': 'right',
                            'width' : '8%'
                        },
                        {
                            'if': {'column_id': 'Gene Names (primary)'},
                            'textAlign': 'left',
                            'width' : '15%'
                        },
                        {
                            'if': {'column_id': 'Organism'},
                            'textAlign': 'left',
                            'width' : '17%'
                        },
                        {
                            'if': {'column_id': 'LeidenCluster_circle'},
                            'textAlign': 'left',
                            'width' : '7%'
                        },
                        {
                            'if': {'column_id': 'Annotation_circle'},
                            'width' : '6%'
                        },
                        {
                            'if': {'column_id': 'Taxon_circle'},
                            'width' : '10%'
                        },
                        {
                            'if': {'column_id': 'Length'},
                            'width' : '6%'
                        }
                    ],
                    tooltip_data = tooltip_data,
                    tooltip_duration=None,
                    filter_action = "native",
                    filter_options = {'case': 'insensitive'},
                    sort_action = "native",

                )]
        )
    ])

    @callback(
        Output('data-table', 'data'),
        Output('data-table', 'tooltip_data'),
        Output('search-bar', 'value'),
        Output('scatter-plot', 'selectedData'),
        Input('search-button', 'n_clicks'),
        Input('reset-button', 'n_clicks'),
        Input('search-bar', 'n_submit'),
        Input('scatter-plot', 'clickData'),
        Input('scatter-plot', 'selectedData'),
        State('search-bar', 'value')
    )
    def search_table(input_button, reset_button, search_bar, clickData, selectedData, value):
        triggered_id = ctx.triggered_id

        if triggered_id == 'search-button' or triggered_id == 'search-bar':
            if value == '':
                return scatter_records, tooltip_data, value, None
            else:
                searched = search_df(scatter_data, value)
                return searched.to_dict('records'), tooltip_maker(searched), value, None
        elif triggered_id == 'reset-button':
            return scatter_records, tooltip_data, '', None
        elif triggered_id == 'scatter-plot':
            if selectedData is not None:
                protids = [i['customdata'][0] for i in selectedData['points'] if 'customdata' in i]
                selected_data = scatter_data[scatter_data['protid'].isin(protids)]

                return selected_data.to_dict('records'), tooltip_maker(selected_data),  '', None
            elif selectedData is None:

                if clickData is not None:
                    protid = clickData['points'][0]['customdata'][0]
                    clicked_data = scatter_data[scatter_data['protid'] == protid]

                    return clicked_data.to_dict('records'), tooltip_maker(clicked_data), '', None
                
                return scatter_records, tooltip_data, '', None
        else:
            return scatter_records, tooltip_data, '', None
        
    @callback(
        Output("scatter-plot", "figure"),
        Input("space-toggle-switch", "value")
    )
    def toggle_space(switch):
        if switch:
            return scatter_umap
        else:
            return scatter_tsne
        
    @callback(
        Output("download-tsv-component", "data"),
        Input("download-button", "n_clicks"),
        State("data-table", "data"),
        prevent_initial_call = True
    )
    def func(n_clicks, data):

        output_data = pd.DataFrame(data)
        output_data = output_data[[col for col in output_data.columns if '_circle' not in col and 'protid_link' not in col]]

        return dcc.send_data_frame(output_data.to_csv, f"{analysis_name}_data.csv")
        
    @callback(
        Output('info-text', 'children'),
        Input('scatter-plot', 'hoverData'),
        Input('tutorial-button', 'n_clicks')
        )
    def display_hover_data(hoverData, tutorial):
        triggered_id = ctx.triggered_id

        color_dict_dict = {
            'Cluster': dict(LC_conditional_color_dict),
            'Score': dict(Annotation_conditional_color_dict),
            'Taxon': taxon_color_dict,
            'Source': source_color_dict,
        }

        tutorial_content = """
            <br>

            <hr/>

            ## **Exploration**

            1. Hover over data in **Scatter** to show protein information.
            2. Click on a point to show it in the **Database** card.
            3. Use the Box Select or Lasso Select tools in the top right of the **Scatter** card to filter the **Database** card.
            4. Export the current contents fo the **Database** card as a CSV usign the **Download** button.

            <hr/>

            ## **Feature Information**
            - **Cluster:** these clusters are identified using the <a href="https://www.cwts.nl/blog?article=n-r2u2a4" target="_blank" children="Leiden algorithm" />.
            - **Score:** the <a href="https://www.uniprot.org/help/annotation_score" target="_blank" children="Uniprot Annotation Score" /> of the protein.
            - **Taxon:** the most-specific taxonomic group as listed in the plot. For example, a human protein be categorized as "Mammalia" because it is more specific than "Vertebrata".
            - **Source:** whether the protein was a hit from BLAST, Foldseek, or both.

        """

        if triggered_id == 'tutorial-button':
            return tutorial_content

        help_text = """
        <br>

        <hr/>

        ## **Welcome to ProteinCartography's Dashboard!**

        Hover over a point to see protein information, or click **Tutorial** above to learn more.
        
        """

        if hoverData is None:
            return help_text
        
        customdata = hoverData['points'][0]['customdata']

        output_dict = {
            'protid': customdata[0],
            'Protein name': customdata[1],
            'Gene name': customdata[2],
            'Organism': customdata[3],
            'Cluster': customdata[4],
            'Score': customdata[5],
            'Taxon': assign_taxon(customdata[6], taxon_order, hierarchical = True),
            'Length': str(customdata[7]) + 'aa',
            'Source': customdata[8]
        }
        
        entries_list = []
        header = ''
        for key,value in output_dict.items():
            if key == 'protid':
                if value in uniprot_set:
                    header = f'<a href="https://www.uniprot.org/uniprotkb/{value}" children="{value}" style="font-size: 20px;" target="_blank"/><hr>'
                else:
                    header = f'<span style="font-size: 20px;" children="{value}" /><hr>'
            elif key in color_dict_dict:
                bg_color = color_dict_dict[key][value]
                bg_color_rgb = mcolors.to_rgb(bg_color)
                color_lightness = colorsys.rgb_to_hls(bg_color_rgb[0], bg_color_rgb[1], bg_color_rgb[2])[1]
                text_color = '#FFFFFF' if color_lightness < 0.8 else '#000000'
                entries_list.append(f'<p><b>{key}:</b> <span style="background-color:{bg_color}; color: {text_color}; padding: 2px; border-radius: 2px;" children="{value}" /></p>')
            else:
                entries_list.append(f'<p><b>{key}:</b> {value}</p>')

        if len(customdata) > 9:
            additional_info = customdata[9:]
            additional_plotting_rules = [i.replace('_v_', ' vs. ').replace('.hit', ' hit') for i in list(plotting_rules.keys())[9:]]
            additional_info = dict(zip(additional_plotting_rules, additional_info))
            additional_md = '<hr>' + '<br>'.join([f'<b>{key}:</b> {value}' for key, value in additional_info.items()])
        else:
            additional_md = ''

        return header + ''.join(entries_list) + additional_md
    
    @callback(
        Output("dashbio-mol3dviewer", "modelData"),
        Output("dashbio-mol3dviewer", "styles"),
        Output("structure-protein-1-input", "value"),
        Input("structure-display-button", "n_clicks"),
        Input("scatter-plot", "clickData"),
        State("structure-protein-1-input", "value"),
        prevent_initial_call = True
    )
    def func(n_clicks, clickData, protein1):
        triggered_id = ctx.triggered_id
        
        if triggered_id == 'scatter-plot':
            if clickData is not None:
                customdata = clickData['points'][0]['customdata'][0]
                protein1 = customdata + '.pdb'

        protein1_path = os.path.join(structures_path, protein1)
        parser1 = PdbParser(protein1_path)

        data1 = parser1.mol3d_data()
        styles1 = create_mol3d_style(
            data1['atoms'], visualization_type='cartoon',
        )
        
        confidences = extract_residue_confidence(protein1_path)
        colors = assign_residue_colors(confidences)

        for i, atom in enumerate(styles1):
            atom['color'] = colors[i]

        return data1, styles1, protein1

    app.run(debug=True)

def main():
    args = parse_args()
    input_path = args.input_path
    analysis_name = args.analysis_name
    keyids = args.keyids
    taxon_focus = args.taxon_focus

    data_dict = make_data_dict(
        input_path = input_path,
        analysis_name = analysis_name, taxon_focus = taxon_focus, keyids = keyids,
    )

    make_dash_app(data_dict)

# Run the app when called from the command line
if __name__ == '__main__':
    main()
