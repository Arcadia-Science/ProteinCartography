#!/usr/bin/env python
import argparse
import textwrap
from typing import Optional, Union

import arcadia_pycolor as apc
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
from color_utils import (
    ANNOTATION_SCORE_COLOR_DICT,
    BAC_COLOR_DICT,
    EUK_COLOR_DICT,
    PDB_ORIGIN_COLOR_DICT,
    PLDDT_CMAP,
    SOURCE_COLOR_DICT,
    arcadia_cividis_r,
    arcadia_magma,
    arcadia_poppies_r,
    arcadia_viridis,
)

# only import these functions when using import *
__all__ = [
    "apply_coordinates",
    "assign_taxon",
    "rescale_list",
    "generate_plotting_rules",
    "plot_interactive",
]

COLORBAR_DEFAULT_DICT = {
    "x": 0.5,
    "y": 0,
    "orientation": "h",
    "xanchor": "center",
    "yanchor": "top",
    "ticklabelposition": "inside top",
    "title_font_size": 14,
    "title_side": "bottom",
    "len": 0.7,
    "thickness": 20,
}

###############
## Functions ##
###############


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dimensions",
        required=True,
        help="path to dimensionality reduction file",
    )
    parser.add_argument("-f", "--features", required=True, help="path to aggregated features file")
    parser.add_argument("-o", "--output", required=True, help="path of destination HTML file")
    parser.add_argument(
        "-t",
        "--dimensions-type",
        help="dimensions data type (pca, tsne, umap, etc.)",
    )
    parser.add_argument(
        "-k",
        "--keyids",
        nargs="*",
        help="keyids for plotting. Usually input protids.",
    )
    parser.add_argument(
        "-x",
        "--taxon-focus",
        default="euk",
        help=(
            "Coloring scheme/ taxonomic groups for broad taxon plot.\n"
            "'euk'(aryote) is default; 'bac'(teria) is another option."
        ),
    )
    parser.add_argument("-X", "--plot-width", default="700", help="width of resulting plot.")
    parser.add_argument("-Y", "--plot-height", default="750", help="width of resulting plot.")
    args = parser.parse_args()

    return args


def rescale_list(values: list[Union[int, float]], new_min: float, new_max: float):
    """
    Rescales a list of values to a new range.

    Args:
        values (list): list of values to rescale.
        new_min (float): minimum value of new range to rescale to.
        new_max (float): maximum value of new range to rescale to.
    """
    # Find the original minimum and maximum values
    original_min = min(values)
    original_max = max(values)

    # Calculate the range of the original values
    original_range = original_max - original_min

    # Calculate the range of the new values
    new_range = new_max - new_min

    # Rescale each value within the new range
    rescaled_values = []
    for value in values:
        # Calculate the relative position of the value within the original range
        relative_position = (value - original_min) / original_range

        # Rescale the relative position to the new range
        rescaled_value = (relative_position * new_range) + new_min

        # Add the rescaled value to the list
        rescaled_values.append(rescaled_value)

    return rescaled_values


def apply_coordinates(
    dimensions_file: str,
    features_file: str,
    saveprefix: Optional[str] = None,
    dimtype: Optional[str] = None,
    save: bool = False,
    prep_step: bool = False,
):
    """
    Adds coordinate information for data from the dimensions file to the aggregated features file.

    Args:
        dimensions_file (str): path to dimensions file.
        features_file (str): path to aggregated features file.
        saveprefix (str): prefix for output file.
        dimtype (str): type of dimensions file (e.g. pca, tsne, etc.).
            If not provided, is inferred from dimensions file.
        save (bool): whether or not to save the file. Defaults to False.
        prep_step (bool): if True, returns path of the final output file. Defaults to False.
    """

    # Load the dimensions data
    reduced_dim_df = pd.read_csv(dimensions_file, sep="\t")

    # If dimtype not provided, infer based on second column of dataframe
    if dimtype is None:
        outdimtype = "".join([i for i in reduced_dim_df.columns[1] if not i.isdigit()])
    else:
        outdimtype = dimtype

    # Check to make sure there's a tsv there
    if ".tsv" not in features_file:
        raise Exception(f'{features_file} does not end in ".tsv" as expected.')

    # Load the features data
    agg_features_df = pd.read_csv(features_file, sep="\t")

    # Merge the two dataframes
    plot_features_df = reduced_dim_df.merge(agg_features_df, on="protid")

    # Generate filepath based on dimensions type and prefix
    if saveprefix is not None:
        savefile = "_".join([saveprefix, outdimtype + ".tsv"])
    else:
        savefile = features_file.replace(".tsv", "_" + outdimtype + ".tsv")

    # Save if needed
    if save:
        plot_features_df.to_csv(savefile, sep="\t", index=None)

    # If prep_step, return the savefile path
    # Otherwise, return the results
    if prep_step:
        return savefile
    else:
        return plot_features_df


def assign_taxon(taxon_list: list, rank_list: list, hierarchical=False, sep=","):
    """
    Takes a list of taxa and assigns it a category based on a ranked-order list.

    The ranked-order list should be ordered from most-specific
    to least-specific taxonomic groupings.
    The function iterates through the list of ranks and checks if each is within the taxon list.
    If hierarchical = True, returns only the first element.
    Otherwise, concatenates the hits in order.

    Args:
        taxon_list (list): list of taxa for the object in question.
        rank_list (list): ranked list to search through.
        hierarchical (bool): whether to return the first element found in the rank_list
            or a concatentation of all hits.
        sep (str): separator character when returning a concatenation of hits.
    Returns:
        a string either containing the first-rank hit or the string concatenation of all hits.
    """
    # check for every element of the rank list whether it's in the taxon list
    output = [item for item in rank_list if item in taxon_list]

    # if output is empty, fill it in with something
    if output == []:
        output.append("Other")

    # if hierarchical, return the lowest rank item
    if hierarchical:
        return output[0]
    # else, return a joined string
    else:
        return sep.join(output)


def generate_plotting_rules(
    taxon_focus: str,
    keyids: Optional[list] = None,
    wordwrap: bool = True,
    num_decimals=4,
) -> dict:
    """
    Generates plotting rules dictionary for use in plot_interactive.

    Args:
        taxon_focus (str): taxonomic focus for broad taxon plot, either 'euk' or 'bac'.
        keyids (list): list of keyids for plotting.
        version (str): version of the aggregated features file.
        wordwrap (bool): whether or not to wordwrap long hovertext fields.
    """

    # if the taxonomic focus is eukaryote, use these groupings and colors
    if taxon_focus == "euk":
        taxon_color_dict = EUK_COLOR_DICT

    # else if the taxonomic focus is bacteria, use these groupings and colors
    elif taxon_focus == "bac":
        taxon_color_dict = BAC_COLOR_DICT

    # use this dictionary as default plotting rules
    # the order of elements in this dictionary determines their order
    # both in the plotting dropdown list and the hovertext
    # long text (>40 characters) in some hovertext fields is now wrapped
    plotting_rules = {
        "Protein names": {
            "type": "hovertext",  # this data only shows up in hovertext
            "fillna": "",
            "textlabel": "Protein name",
        },
        "Gene Names (primary)": {
            "type": "hovertext",  # this data only shows up in hovertext
            "fillna": "",
            "textlabel": "Gene name",
        },
        "Organism": {
            "type": "hovertext",  # this data only shows up in hovertext
            "fillna": "",
            "textlabel": "Organism",
        },
        "LeidenCluster": {
            "type": "categorical",
            "fillna": "None",
            "apply": lambda x: str(
                x
            ),  # Leiden cluster is often read as int; this forces it to be string
            "color_order": apc.Palettes["arcadia:AccentAllOrdered"].colors,
            "textlabel": "Leiden Cluster",
        },
        "Annotation": {
            "type": "categorical",
            "fillna": 0,
            "apply": lambda x: str(
                int(x)
            ),  # Annotation score is parsed as float but we want it to be string
            "color_dict": ANNOTATION_SCORE_COLOR_DICT,
            "textlabel": "Annotation Score",
        },
        "Lineage": {
            "type": "taxonomic",
            "fillna": "[]",
            "apply": lambda x: eval(
                x
            ),  # This converts the taxonomic groupings from a string-ified list to a real list
            "taxon_order": taxon_color_dict.keys(),
            "color_order": taxon_color_dict.values(),
            "textlabel": "Broad Taxon",
            "skip_hover": True,
        },
        "Length": {
            "type": "continuous",
            "fillna": 0,
            "textlabel": "Length",
            "color_scale": arcadia_cividis_r,
        },
        "source.method": {
            "type": "categorical",
            "fillna": "None",
            "color_dict": SOURCE_COLOR_DICT,
            "textlabel": "Source",
        },
        "pdb_origin": {
            "type": "categorical",
            "fillna": "None",
            "color_dict": PDB_ORIGIN_COLOR_DICT,
            "textlabel": "PDB Origin",
        },
        "pdb_confidence": {
            "type": "continuous",
            "fillna": -1,
            "textlabel": "Mean pLDDT",
            "color_scale": PLDDT_CMAP,
            "cmin": 0,
            "cmax": 100,
            "apply": lambda x: round(x, num_decimals),
        },
    }

    # Add additional plotting rules for each keyid (usually the input protein)
    if keyids is not None:
        for keyid in keyids:
            # Add a plot for tmscore to each input protein
            plotting_rules[f"TMscore_v_{keyid}"] = {
                "type": "continuous",
                "apply": lambda x: round(x, num_decimals),
                "fillna": -0.01,
                "textlabel": f"TMscore vs. {keyid}",
                "color_scale": arcadia_viridis,
                "cmin": 0,
                "cmax": 1,
            }
            plotting_rules[f"fident_v_{keyid}"] = {
                "type": "continuous",
                "apply": lambda x: round(x, num_decimals),
                "fillna": -0.01,
                "textlabel": f"Fraction seq identity vs. {keyid}",
                "color_scale": arcadia_magma,
                "cmin": 0,
                "cmax": 1,
            }
            plotting_rules[f"concordance_v_{keyid}"] = {
                "type": "continuous",
                "apply": lambda x: round(x, num_decimals),
                "fillna": -1.01,
                "textlabel": f"Concordance vs. {keyid}",
                "color_scale": arcadia_poppies_r,
                "cmin": -1,
                "cmax": 1,
            }
            # Add hovertext for whether or not a given protein was a hit via blast or foldseek
            # to the input protein
            plotting_rules[f"{keyid}.hit"] = {
                "type": "hovertext",
                "apply": lambda x: True if x == 1 else False,
                "fillna": "None",
                "textlabel": f"Blast/Foldseek Hit to {keyid}",
            }

    # wordwrap long columns in hovertext
    if wordwrap:
        for rule in plotting_rules:
            if rule in ["Protein names", "Gene Names (primary)", "Organism"]:
                plotting_rules[rule]["apply"] = (
                    lambda x: "</br>".join(textwrap.wrap(x, width=40)),
                )

    return plotting_rules


def preprocess_dataframe(
    plotting_rules: dict,
    df: pd.DataFrame,
):
    """
    Preprocesses a dataframe based on the plotting rules.
    Fills NAs and applies functions as needed.

    Args:
        plotting_rules (dict): plotting rules dictionary.
        df (pd.DataFrame): dataframe to preprocess.
    """
    for col in plotting_rules:
        if col not in df.columns:
            continue

        # if the plotting rule 'fillna' is present, fills NAs with that value
        if "fillna" in plotting_rules[col].keys():
            df[col] = df[col].fillna(plotting_rules[col]["fillna"])

        # if the plotting rule 'apply' is present, applies that function
        if "apply" in plotting_rules[col].keys():
            df[col] = df[col].apply(plotting_rules[col]["apply"])


def generate_hovertemplate(
    plotting_rules: dict,
    df: pd.DataFrame,
    custom_data: list,
):
    """
    Generates a hovertemplate based on the plotting rules.
    Plotly hovertemplates are a text string that can be formatted with custom data.

    Args:
        plotting_rules (dict): plotting rules dictionary.
        df (pd.DataFrame): dataframe to use for generating hovertemplate.
    """
    # generate a hover template text string
    # custom data in plotly is indexed by order
    hovertemplate_generator = ["<b>%{customdata[0]}</b></br>––––––––––––"]
    # pass in the custom data columns
    # this sets protid as customdata[0]
    # and the rest of the plotting rules in order

    # iterate through plotting rules, applying relevant rules
    for col in plotting_rules:
        if col not in df.columns:
            continue

        # if you want to skip an element from being added to hover text
        if "skip_hover" in plotting_rules[col].keys():
            continue

        # sets up hover text based on 'textlabel' attributes
        if "textlabel" in plotting_rules[col].keys():
            hovertext_label = plotting_rules[col]["textlabel"]
        else:
            hovertext_label = col

        # generates hoverlabel custom text string for that column
        # This doesn't work properly if the column is NA for that value.
        col_index = custom_data.index(col)

        hovertemplate_item = "<b>" + hovertext_label + "</b>: %{customdata[" + str(col_index) + "]}"
        hovertemplate_generator.append(hovertemplate_item)

    hovertemplate = "<br>".join(hovertemplate_generator)

    return hovertemplate


def generate_scatterplot(
    plotting_rule: dict,
    col: str,
    df: pd.DataFrame,
    custom_data: list,
    dim1: str,
    dim2: str,
):
    """
    Generates a scatterplot based on the plotting rules.

    Args:
        plotting_rule (dict): plotting rule dictionary for the specific column.
        col (str): column name to plot based on rules.
        df (pd.DataFrame): dataframe to use for generating hovertemplate.
        custom_data (list): list of custom data columns that were included in the hovertemplate.
        dim1 (str): name of first dimension, the x value.
        dim2 (str): name of second dimension, the y value.
    """
    # Plotting rules for categorical plots
    if plotting_rule["type"] == "categorical":
        color_keys = np.sort(df[col].unique())

        # if a color dict is specified, use that color dict
        if "color_dict" in plotting_rule.keys():
            colors_dict = plotting_rule["color_dict"]

            # make sure there's at least one entry for every possible value for that column
            # in the future, should automatically add extra colors.
            if not all([entry in colors_dict.keys() for entry in color_keys]):
                raise Exception("color_dict is missing some entries.")

        # if a color dict is not provided, but a color order is, use that to make a color dict
        elif "color_order" in plotting_rule.keys():
            color_order = list(plotting_rule["color_order"])

            # TODO: consider whether we should support this use case
            # TODO: maybe we serialize a palette and save it to a file?

            # make sure there's enough colors to make a dict with
            # if not, extend available colors
            if len(color_order) < len(color_keys):
                more_colors = apc.extend_colors(color_order, len(color_keys))

                # extend color order with new colors
                color_order = more_colors

            # make color dict
            colors_dict = dict(zip(color_keys, color_order))

        # set remaining variables to pass to the plotting function
        category_orders = {col: color_keys}
        color_scale = None

    # Plotting rules for continuous plots
    elif plotting_rule["type"] == "continuous":
        if "color_scale" in plotting_rule.keys():
            color_scale = plotting_rule["color_scale"]
        else:
            color_scale = "viridis"

        # set remaining variables to pass to the plotting function
        category_orders = None
        colors_dict = None

    # Plotting rules for taxonomic plots
    elif plotting_rule["type"] == "taxonomic":
        # Check to make sure there's a taxon order
        if "taxon_order" in plotting_rule.keys():
            taxon_order = plotting_rule["taxon_order"]
        else:
            raise Exception('Please provide a "taxon_order" list in the plotting rules.')

        # Set color keys to taxon order
        color_keys = taxon_order

        # if a color dict is specified, use that color dict
        if "color_dict" in plotting_rule.keys():
            colors_dict = plotting_rule["color_dict"]

            # make sure there's at least one entry for every possible value for that column
            # in the future, should automatically add extra colors.
            if not all([entry in colors_dict.keys() for entry in color_keys]):
                raise Exception("color_dict is missing some entries.")

        # if a color dict is not provided, but a color order is, use that to make a color dict
        elif "color_order" in plotting_rule.keys():
            color_order = plotting_rule["color_order"]

            # make sure there's enough colors to make a dict with
            # if not, extend available colors
            if len(color_order) < len(color_keys):
                # extend color order with new colors
                color_order = apc.extend_colors(color_order, len(color_keys))

            # make colors dict
            colors_dict = dict(zip(color_keys, color_order))

            # add entry for 'Other'
            colors_dict["Other"] = apc.All["arcadia:brightgrey"]

        # generate a taxonomic category column using assign_taxon lambda function
        col_taxonomic = col + "_taxonomic"
        df[col_taxonomic] = df[col].apply(
            lambda x, taxon_order=taxon_order: assign_taxon(x, taxon_order, hierarchical=True)
        )
        col = col_taxonomic
        category_orders = {col: color_keys}
        color_scale = None

    # generate a plot using the above parameters
    plot = px.scatter(
        df,
        dim1,
        dim2,
        color=col,
        hover_name="protid",
        category_orders=category_orders,
        color_discrete_map=colors_dict,
        color_continuous_scale=color_scale,
        custom_data=custom_data,
    )
    return plot


def regenerate_scatterplot(
    plots: dict,
    plotting_rules: dict,
    df: pd.DataFrame,
    scatter_counter: dict,
    marker_opacity: float = 0.8,
    marker_size: int = 4,
):
    """
    Regenerates a scatterplot based on the plotting rules.

    Because we want to have a dropdown menu to toggle visibility of different plots,
    where each plot has variable numbers of traces,
    and some plots are categorical while others are continuous,
    and we want to preserve the click/hover functionality of the plot,
    we need to hack around Plotly's limitations.

    Building initial plots as individual px.scatter objects, then transferring them to a new figure
    is a way we can achieve this.

    This allows us to have custom color scales, cmin, cmax,
    and control visibility for groups of traces using a dropdown menu.

    TODO: There may be another way to solve this, but this is the solution that works for now.

    Args:
        plots (dict): dictionary of plotly scatter objects to transfer to a new figure.
        plotting_rules (dict): full plotting rules dictionary.
        df (pd.DataFrame): dataframe to use for generating hovertemplate.
        marker_opacity (float): opacity of markers  to be passed to go.Scatter or go.Scattergl.
        marker_size (int): size of markers to be passed to go.Scatter or go.Scattergl.
    """
    # create a blank figure
    fig = go.Figure()

    # iterate through the plots to get their objects
    for j, (col, plot) in enumerate(plots.items()):
        # hide all plots except the first one by default
        if j == 0:
            vis = True
        else:
            vis = False

        # if there are more than 50 different categories, hide the legend
        # this keeps the plot's aspect ratio from changing too much
        if scatter_counter[col] > 50:
            showlegend = False
        else:
            showlegend = True

        # transfer each individual plot's traces to the new empty figure
        for scatter in plot.data:
            # detect whether a trace is a Scattergl or Scatter object
            obj_method = go.Scatter
            if isinstance(scatter, plotly.graph_objs._scattergl.Scattergl):
                obj_method = go.Scattergl

            # collect color bar style parameters
            colorbar_dict = COLORBAR_DEFAULT_DICT
            colorbar_dict["title"] = col

            # Setting `cmin` and `cmax` when using px.scatter requires using `update_layout`,
            # which may not transfer properly to the new figure.
            # Also, because the min and max values for each scatter plot are different,
            # we can't set it globally.
            # So, we need to set it manually upon adding traces to the new figure.
            # TODO: There may be another way to solve this,
            # but this is the solution that works for now.
            if plotting_rules[col]["type"] == "continuous":
                # rescale the colorbar to user-defined values if needed
                # this is needed when (e.g. for pLDDT scores) the maximum value computed
                # is not the maximum value possible
                if "cmax" in plotting_rules[col]:
                    cmax = plotting_rules[col]["cmax"]
                else:
                    cmax = df[col].max()

                if "cmin" in plotting_rules[col]:
                    cmin = plotting_rules[col]["cmin"]
                else:
                    cmin = df[col].min()

                if "color_scale" in plotting_rules[col]:
                    new_color_scale = plotting_rules[col]["color_scale"]
                else:
                    new_color_scale = "viridis"

                if "fillna" in plotting_rules[col] and plotting_rules[col]["fillna"] < cmin:
                    fillna_value = plotting_rules[col]["fillna"]
                    fillna_fraction = -1 * (cmin - fillna_value) / (cmax - fillna_value)

                    input_values = [fillna_fraction] + [i[0] for i in new_color_scale]
                    original_colors = [i[1] for i in new_color_scale]
                    new_values = rescale_list(input_values, 0, 1)
                    new_values_discretized = new_values[0:2] + new_values[1:]

                    if "na_color" in plotting_rules[col]:
                        na_color = plotting_rules[col]["na_color"]
                    else:
                        na_color = apc.All["arcadia:brightgrey"]

                    new_colors = [na_color] * 2 + original_colors

                    new_color_scale_collector = [
                        [new_values_discretized[i], new_colors[i]]
                        for i in np.arange(len(new_values_discretized))
                    ]
                    new_color_scale = new_color_scale_collector

                    cmin = plotting_rules[col]["fillna"]

                fig.add_trace(
                    obj_method(
                        scatter,
                        visible=vis,
                        marker=dict(
                            color=scatter.marker.color,
                            colorbar=colorbar_dict,
                            colorscale=new_color_scale,
                            opacity=marker_opacity,
                            size=marker_size,
                            cmax=cmax,
                            cmin=cmin,
                        ),
                    )
                )
            else:
                fig.add_trace(obj_method(scatter, visible=vis, showlegend=showlegend))
    return fig


def generate_dropdown(scatter_counter, plotting_rules, plots):
    # define a helper function to determine the visibility toggle list
    # this sets, for every single trace (each category of each plot gets a trace), its visibility
    def visibility_list(col):
        entries = []

        for fig in scatter_counter.keys():
            if fig == col:
                entries.extend([True] * scatter_counter[fig])
            else:
                entries.extend([False] * scatter_counter[fig])

        return entries

    # collector for buttons
    buttons = []

    # for each plot, create a button
    for col in plots.keys():
        # get the textlabel for that plot if it's provided
        # otherwise, default to column name
        if "textlabel" in plotting_rules[col].keys():
            button_label = plotting_rules[col]["textlabel"]
        else:
            button_label = col

        # create a button that toggles visibility based on which data is currently selected
        button_item = dict(
            args=["visible", visibility_list(col) + [True]],
            label=button_label,
            method="restyle",
        )
        # add this button to the selector
        buttons.append(button_item)

    # collector for dropdown menu
    dropdown_menu = dict(
        buttons=list(buttons),
        showactive=True,
        x=0.08,
        xanchor="left",
        y=1.02,
        yanchor="bottom",
        font_size=14,
        bgcolor="white",
    )

    return dropdown_menu


def plot_interactive(
    coordinates_file: str,
    plotting_rules: dict,
    marker_size: int = 4,
    marker_opacity: float = 0.8,
    output_file: Optional[str] = None,
    keyids: Optional[list] = None,
    show: bool = False,
    plot_width: int = 700,
    plot_height: int = 750,
    plot_bgcolor=apc.All["arcadia:paper"],
    paper_bgcolor="rgba(0,0,0,0)",
    hide_hover: bool = False,
):
    """
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
        - If 'hovertext', won't plot the data as a dropdown menu option,
            but will include the data in the hover-over text box.
    - **'fillna':**
        - A value to fill rows that are `np.nan`.
        - For categorical plots, usually empty string `''`.
        - For continuous plots, usually `0`.
    - **'apply':**
        - A function that will be applied to every element of the feature before plotting.
        - This can be used to convert a numerical value into a categorical value
            with a lambda function.
        - For example, `lambda x: str(x)`
    - **'skip_hover':**
        - Boolean, whether to include this data in the hovertext.
    - **'textlabel':**
        - A string that replaces the column name on buttons and in hover text.
        - Useful if the column name is too ugly.

    ### For 'categorical' and 'taxonomic' plots ###
    - **'color_order':**
        - A list of HEX code colors used for coloring categorical variables.
        - If there aren't enough colors for unique values of the data,
            will generate up to 3x more colors of varying darkness.
    - **'color_dict':**
        - A dictionary of key: value pairs for coloring categorical variables.
        - The key is the name of the category and the value is a HEX code color.

    ### For 'taxonomic' plots ###
    - **'taxon_order':**
        - Exclusively used for 'taxonomic' style plots.
            A list of ranked-order taxa for categorization.
        - Should start with more-specific taxa and expand to less-specific taxa.

    Args:
        coordinates_file (str): path to coordinates file (an aggregated features file
            with columns for the coordinates in PC/TSNE/UMAP space).
        plotting_rules (dict): a JSON-like dictionary containing plotting rules
            for relevant columns of the data.
        marker_size (int): size of markers.
        marker_opacity (float): opacity of markers (value should between 0 and 1).
        output_file (str): path of destination file.
        keyids (list): list of key protids to assign a star-shaped marker. Usually input proteins.
        show (bool): whether or not to show the plot.
        hide_hover (bool): whether to override the hover text and just show the protid.
    Returns:
        if show = False, returns the plotly.graphobjects object of the plot.
    """

    # make a copy of the starting data
    df = pd.read_csv(coordinates_file, sep="\t")

    dim1 = df.columns[1]
    dim2 = df.columns[2]

    # preprocesses the dataframe based on the plotting rules
    preprocess_dataframe(plotting_rules, df)

    valid_plotting_rules = [col for col in plotting_rules.keys() if col in df.columns]
    custom_data = ["protid"] + list(valid_plotting_rules)

    # generates a full hovertemplate string from hovertemplate_generator list
    hovertemplate = generate_hovertemplate(plotting_rules, df, custom_data)

    # if hide_hover is true, override the hovertemplate
    if hide_hover:
        hovertemplate = "<b>%{customdata[0]}</b>"

    # Collector dictionary for making plots
    plots = {}

    # Iterate through the plotting rules for each plotted datatype
    # We make each plot for each datatype individually;
    # Then, we transfer the points from the existing plots to a new, single plot
    # This way we can get the toggle system working
    for col in plotting_rules.keys():
        # ignore invalid rules
        if col not in df.columns:
            continue

        plotting_rule = plotting_rules[col]
        if plotting_rule["type"] == "hovertext":
            continue

        # Generate a plot based on rules
        plots[col] = generate_scatterplot(
            plotting_rules[col],
            col,
            df,
            custom_data,
            dim1,
            dim2,
        )
        # Update the plot with the hovertemplate
        plots[col].update_traces(
            marker=dict(size=marker_size, opacity=marker_opacity),
            hovertemplate=hovertemplate,
        )

    # keep track of how many traces there are per plot
    # for example, a categorical trace would have one trace per category
    # whereas a continuous trace would have one trace total
    # this is used to toggle visibility of traces using a dropdown later
    scatter_counter = {}

    # iterate through the plots to get their objects
    for col, plot in plots.items():
        # counts the number of different scatter traces within the original plot
        # this is needed to determine the number of true or false values
        # for the dropdown toggle visibility field
        scatter_counter[col] = len(plot.data)

    # create new empty figure to move every original figure onto
    # this is a workaround to allow colored legends to be preserved
    # and automatically switched using the dropdown
    # the alternative approach forces you to show all legend items at the same time
    # for all colorations, which is messy
    fig = regenerate_scatterplot(
        plots,
        plotting_rules,
        df,
        scatter_counter,
        marker_opacity=marker_opacity,
        marker_size=marker_size,
    )

    # generate the dropdown menu with dynamic toggle visibility per trace for each plotting rule
    dropdown_menu = generate_dropdown(
        scatter_counter,
        plotting_rules,
        plots,
    )

    # if there are any keyids provided, generate an additional plot and add a toggle button
    if keyids is not None:
        # get the positions of the key ids
        keypoints = df[df["protid"].isin(keyids)]

        # add those key points as a new scatter with a star-diamond marker
        fig.add_trace(
            go.Scatter(
                mode="markers",
                x=keypoints[dim1],
                y=keypoints[dim2],
                marker=dict(color="rgba(0,0,0,0.9)", size=12, symbol="star-diamond"),
                showlegend=False,
                hoverinfo="skip",
                hoverlabel=None,
                name="keyids",
            )
        )

        # create a separate button to toggle keyid trace
        # the dictionary within `buttons` allows you to modify just the last trace in the plot
        keyids_button = dict(
            buttons=[
                dict(
                    # [parameter to modify, value to modify it at]
                    args=[{"visible": True}, [len(fig.data) - 1]],
                    args2=[{"visible": False}, [len(fig.data) - 1]],
                    label="Input Proteins",
                )
            ],
            type="buttons",
            showactive=True,
            x=1,
            xanchor="right",
            y=1.02,
            yanchor="bottom",
            font_size=14,
        )
        # add this to the menu buttons
        updatemenu = [dropdown_menu, keyids_button]
    else:
        # if no keyids, don't add this button
        updatemenu = [dropdown_menu]

    # add the menu buttons
    fig.update_layout(updatemenus=updatemenu)

    # set the width, height, and other parameters of the plot
    fig.update_layout(
        width=plot_width,
        height=plot_height,
        annotations=[
            dict(
                text="color",
                x=0,
                xref="paper",
                xanchor="left",
                y=1.03,
                yref="paper",
                yanchor="bottom",
                align="right",
                showarrow=False,
                font_size=14,
            )  # This shows the word "color" next to the dropdown
        ],
        plot_bgcolor=plot_bgcolor,
        paper_bgcolor=paper_bgcolor,
    )

    # fix the scale of the x and y axes across plots
    xmin = df[dim1].min()
    xmax = df[dim1].max()
    xwiggle = 0.07 * (xmax - xmin)

    ymin = df[dim2].min()
    ymax = df[dim2].max()
    ywiggle = 0.07 * (ymax - ymin)

    fig.update_layout(margin=dict(l=10, r=10, t=75, b=100))
    # make x and y axes have a fixed scale
    fig.update_yaxes(
        range=(ymin - ywiggle, ymax + ywiggle),
        showticklabels=False,
        showspikes=True,
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="#EAEAEA",
        mirror=True,
    )
    fig.update_xaxes(
        range=(xmin - xwiggle, xmax + xwiggle),
        showticklabels=False,
        showspikes=True,
        showgrid=False,
        zeroline=False,
        showline=True,
        linecolor="#EAEAEA",
        mirror=True,
    )

    # move the legend for categorical and color bar features to the bottom
    fig.update_layout(
        legend=dict(
            orientation="h",
            yanchor="top",
            y=0,
            xanchor="left",
            x=0,
            font=dict(size=12),
            tracegroupgap=1,
        )
    )

    # adjust modebar parametrs
    fig.update_layout(
        modebar={
            "bgcolor": "rgba(0,0,0,0)",
            "color": apc.All["arcadia:denim"],
            "activecolor": apc.All["arcadia:marineblue"],
        }
    )

    # set configurations for scatter plot
    scatter_config = {
        "displayModeBar": True,
        "scrollZoom": True,
        "toImageButtonOptions": {"filename": "scatter", "format": "svg"},
        "modeBarButtonsToRemove": ["zoomIn", "zoomOut"],
    }

    # Use Arial as the default font if possible
    try:
        fig.update_layout(font=dict(family="Arial"))
    except Exception:
        pass

    # save if filename is provided
    if output_file is not None:
        fig.write_html(output_file, config=scatter_config)

    # show if desired
    if show:
        fig.show(config=scatter_config)

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

    # generate coordinates file for the plot
    coordinates_file = apply_coordinates(
        dimensions_file,
        features_file,
        save=True,
        prep_step=True,
        dimtype=dimensions_type,
    )

    # generate plotting rules
    plotting_rules = generate_plotting_rules(taxon_focus, keyids)

    # make the plot
    plot_interactive(
        coordinates_file,
        plotting_rules,
        output_file=output_file,
        keyids=keyids,
        plot_width=plot_width,
        plot_height=plot_height,
    )


# check if called from interpreter
if __name__ == "__main__":
    main()
