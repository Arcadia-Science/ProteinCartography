#!/usr/bin/env python
import argparse
import textwrap

import arcadia_pycolor as apc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.ticker import MaxNLocator
from plot_interactive import extend_colors
from plotly.subplots import make_subplots
from wordcloud import WordCloud

__all__ = [
    "plot_semantic_analysis",
    "semantic_barchart_plotly",
    "wordcloud_image",
    "semantic_multiplot_plotly",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--features-file",
        required=True,
        help="Path to features file for grouping.",
    )
    parser.add_argument("-c", "--agg-column", required=True, help="Column to aggregate groups on.")
    parser.add_argument(
        "-n", "--annot-column", required=True, help="Column of annotations to analyze."
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Path to output file ending in a matplotlib-compatible format "
            "(png, pdf, ps, eps, and svg)."
        ),
    )
    parser.add_argument(
        "-i",
        "--interactive",
        help='Path to output file ending with ".html" for interactive plot.',
    )
    parser.add_argument(
        "-e",
        "--exclude-words",
        nargs="+",
        default=["protein", "None"],
        help="Words to exclude from word cloud.",
    )
    parser.add_argument("-a", "--analysis-name", help="name of analysis for plotting")
    args = parser.parse_args()

    return args


def plot_semantic_analysis(
    features_file: str,
    agg_col: str,
    annot_col: str,
    colors: list,
    exclude_words=("protein", "None"),
    ignore_nan=True,
    top_n=10,
    max_str_len=30,
    n_cols=3,
    analysis_name="",
    savefile=None,
    show=False,
):
    """
    Takes a features file and performs semantic analysis on an annot_col
    for a given group of entries specified in agg_col.
    The annot_col usually contains a gene description for each protein in the features file.
    The agg_col is a column that describes the group for each protein
    (e.g. its cluster number, species, etc).

    Generates a bar chart where the most common annotations (up to top_n) in each group are ranked.
    Also generates a word cloud for the most common words in the annotations for that group.

    Args:
        features_file (str): path of input features file.
        agg_col (str): column in the features file to use for aggregation.
        annot_col (str): column in the features file that contains the annotations.
        colors (list): ordered list of HEX or rgba values to color the groups with.
        exclude_words (list): list of words to ignore when building word cloud.
            Defaults to 'protein' and 'None'.
        ignore_nan (bool): whether to ignore empty annotation cells
            when building the bar chart and word cloud. Defaults to True.
        top_n (int): up to this number of annotations will be displayed in the bar chart.
        max_str_len (int): maximum number of characters to display
            from annotation string in bar chart. Defaults to 30.
        n_cols (int): number of columns of paired bar chart + word cloud plots to show.
            Number of rows is automatically adjusted to fit.
        analysis_name (str): an analysis name to add to the title of the plot. Defaults to ''.
        savefile (str or None): if not None, saves a file to this path.
        show (bool): whether or not to show the plot. Used during interactive sessions.
            Defaults to False.
        output_file (str): path of destination file.
    """
    # read in features file
    features_df = pd.read_csv(features_file, sep="\t")

    # TODO (KC): this is duplicated in count_
    def ignore_function(x):
        if ignore_nan:
            return [i for i in x if i is not np.nan]
        return list(x)

    # group features file by aggregation column and extract aggregated annotation column
    groupedby_agg_df = features_df.groupby(agg_col).agg(ignore_function)[annot_col]

    # determine number of groups
    n_groups = len(groupedby_agg_df)

    used_colors = colors

    if len(used_colors) < n_groups:
        used_colors = used_colors + extend_colors([1] * n_groups, used_colors)

    # set plot row parameters based on number of groups and columns
    n_rows = int(np.ceil(len(groupedby_agg_df) / n_cols))

    # collectors for plot information
    summary_dict = {}
    str_summary_dict = {}
    len_dict = {}
    wc_dict = {}

    # generate summary statistics
    for i, (clu, values) in enumerate(groupedby_agg_df.items()):
        # count the number of occurrences of each exact annotation string
        summary_dict[clu] = pd.DataFrame(pd.value_counts(values))

        # count number of unique annotations per cluster
        len_dict[clu] = len(values)

        # combine all annotations into one long space-separated string,
        # then break into individual words
        annot_word_list = " ".join(list(values)).split(" ")

        # sanitize word list by removing irrelevant words and parentheses
        sanitized_word_list = [
            word.replace("(", "").replace(")", "")
            for word in annot_word_list
            if word not in exclude_words
        ]

        # get value counts per-word
        str_summary = dict(pd.value_counts(sanitized_word_list, normalize=True))

        # save word frequencies to dict
        str_summary_dict[clu] = str_summary

        # generate word cloud based on frequencies
        wc_dict[clu] = WordCloud(
            width=500,
            height=500,
            background_color="white",
            color_func=lambda *args, i=i, **kwargs: used_colors[i],
        ).generate_from_frequencies(str_summary)

    # create figure with correct number of dimensions
    plt.figure(figsize=(n_cols * 6, n_rows * 3))
    plt.suptitle(f"Simple semantic analysis of {analysis_name} {agg_col}", y=1.01, fontsize=18)

    # generate plots
    for i, clu in enumerate(summary_dict.keys()):
        # plot the bar chart
        plt.subplot(n_rows, n_cols * 2, i * 2 + 1)
        top_n_df = summary_dict[clu].head(top_n)

        # shorten annotation strings based on a maximum number of characters
        labels = [
            i[:max_str_len] + "..." if len(i) > max_str_len else i for i in list(top_n_df.index)
        ]
        widths = list(top_n_df["count"])

        plt.barh(
            y=np.arange(len(widths)),
            width=widths,
            tick_label=labels,
            color=used_colors[i],
            alpha=0.8,
        )
        plt.gca().invert_yaxis()
        plt.setp(plt.gca().yaxis.get_majorticklabels(), ha="left", x=0.05)
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.gca().tick_params(axis="y", length=0)
        plt.ylabel(f"{agg_col} {clu}", fontsize=14)
        plt.xlabel("number of annotations")
        plt.title(f"top {top_n} full annotations")

        # plot the word cloud
        plt.subplot(n_rows, n_cols * 2, i * 2 + 2)
        plt.imshow(wc_dict[clu])

        # hide xticks and yticks for word cloud
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.title("proportional word cloud")

    # tighten up the layout after plotting
    plt.tight_layout()

    if savefile is not None:
        plt.savefig(savefile, bbox_inches="tight")

    if show:
        plt.show()


def count_features(
    features_file: str,
    agg_col: str,
    annot_col: str,
    colors: list,
    exclude_words=("protein", "None"),
    ignore_nan=True,
):
    """
    Generate feature counts for a features file.

    Args:
        features_file (str): path to a features file
        agg_col (str): which column to group by
        annot_col (str): which column's annotations to summarize
        colors (list): list of HEX codes to use for the SVG images
        exclude_words (list): words to mask out from annotations
        ignore_nan (bool): whether to ignore NaN annotations
    """
    # read in features file
    features_df = pd.read_csv(features_file, sep="\t")

    # TODO (KC): this is duplicated from plot_semantic_analysis above
    def ignore_function(x):
        if ignore_nan:
            return [i for i in x if i is not np.nan]
        return list(x)

    # group features file by aggregation column and extract aggregated annotation column
    groupedby_agg_df = features_df.groupby(agg_col).agg(ignore_function)[annot_col]

    # determine number of groups
    n_groups = len(groupedby_agg_df)

    used_colors = colors

    if len(used_colors) < n_groups:
        used_colors = used_colors + extend_colors([1] * n_groups, used_colors)

    # collectors for plot information
    annotation_count_dict = {}
    str_annotation_count_dict = {}
    total_annots_dict = {}
    wordclouds_dict = {}

    # generate summary statistics
    for i, (clu, values) in enumerate(groupedby_agg_df.items()):
        # count the number of occurrences of each exact annotation string
        annotation_count_dict[clu] = pd.DataFrame(pd.value_counts(values))

        # count number of unique annotations per cluster
        total_annots_dict[clu] = len(values)

        # combine all annotations into one long space-separated string,
        # then break into individual words
        annot_word_list = " ".join(list(values)).split(" ")

        # sanitize word list by removing irrelevant words and parentheses
        sanitized_word_list = [
            word.replace("(", "").replace(")", "")
            for word in annot_word_list
            if word not in exclude_words
        ]

        # get value counts per-word
        str_summary = dict(pd.value_counts(sanitized_word_list, normalize=True))

        # save word frequencies to dict
        str_annotation_count_dict[clu] = str_summary

        # generate word cloud based on frequencies
        wordclouds_dict[clu] = WordCloud(
            width=500,
            height=500,
            background_color="white",
            color_func=lambda *args, i=i, **kwargs: used_colors[i],
        ).generate_from_frequencies(str_summary)

    results = {
        "annotation_count": annotation_count_dict,
        "str_annotation_count": str_annotation_count_dict,
        "total_annots": total_annots_dict,
        "wordclouds": wordclouds_dict,
    }

    return results


def semantic_barchart_plotly(annotation_count: dict, group: str, color: str, top_n=10):
    """
    Generate a plotly barchart object for the top n annotations

    Args:
        annotation_count (dict): annotation count dictionary output of count_features()
        group (str): which subgroup of agg_col to use for filtering
        color (str): HEX color for that group
        top_n (int): number of annotations to display
    """
    annotation_count_group = annotation_count[group]

    annotation_count_group_filtered = annotation_count_group.head(top_n)

    x = annotation_count_group_filtered["count"]
    text = annotation_count_group_filtered.index
    customdata = ["<br>".join(textwrap.wrap(i, 30)) for i in text]

    bar = go.Bar(
        x=x,
        text=text,
        marker_color=color,
        xaxis="x",
        yaxis="y",
        customdata=customdata,
        hovertemplate="<br>".join(["%{customdata}", "<b>Count:</b> %{x}"]) + "<extra></extra>",
    )
    return bar


def wordcloud_image(wordclouds: dict, group: str, color: str, mode="fig", savefile=None):
    """
    Generate an SVG image or Plotly figure object from a wordcloud object

    Args:
        wordclouds (dict): wordcloud list output of count_features()
        group (str): which subgroup of agg_col to use for filtering
        color (str): HEX color for that group
        mode (str): 'fig' for Plotly object or 'svg'/'png' to save as an image
        savefile (str): path of destination save file when using 'svg' or 'png'
    """

    if mode == "fig":
        wc = wordclouds[group].to_array()
        image = go.Image(z=wc, xaxis="x2", yaxis="y2", hoverinfo="skip")
        return image

    elif mode == "svg" and savefile is not None:
        wordclouds[group].to_svg(savefile)

    elif mode == "png" and savefile is not None:
        wordclouds[group].to_file(savefile)

    return


def semantic_multiplot_plotly(
    count_features_results: dict,
    colors: list,
    n_cols=3,
    savefile=None,
    show=False,
):
    """
    Generate an multiple bar chart/ word cloud chart.

    Args:
        count_features_results (dict): full dictionary output of count_features()
        colors (list): list of HEX colors
        n_cols (int): number of columns to use in plot
        savefile (str): path to destination file
        show (bool): whether to show the result in an interactive session
    """

    n_groups = len(count_features_results["annotation_count"].keys())

    # set plot row parameters based on number of groups and columns
    n_rows = int(np.ceil(n_groups / n_cols))

    fig = make_subplots(
        rows=n_rows,
        cols=n_cols * 2,
        horizontal_spacing=0.02,
        vertical_spacing=0.08,
        subplot_titles=[
            item
            for key in count_features_results["annotation_count"].keys()
            for item in (key, "proportional word cloud")
        ],
    )

    flattened_indices = [
        k
        for lst in [[(j + 1, i + 1) for i in np.arange(n_cols * 2)] for j in np.arange(n_rows)]
        for k in lst
    ]

    used_colors = colors

    if len(used_colors) < n_groups:
        used_colors = used_colors + extend_colors([1] * n_groups, used_colors)

    colors_dict = dict(zip(count_features_results["annotation_count"].keys(), used_colors))

    xaxis_params = {
        "showline": True,
        "linewidth": 1,
        "linecolor": apc.All["arcadia:crow"],
        "title": "number of annotations",
        "title_standoff": 2,
    }

    yaxis_params = {
        "showline": True,
        "linewidth": 1,
        "linecolor": apc.All["arcadia:crow"],
        "autorange": "reversed",
        "showticklabels": False,
    }

    xaxis2_params = {
        "showticklabels": False,
        "showline": True,
        "linewidth": 1,
        "linecolor": apc.All["arcadia:brightgrey"],
        "mirror": True,
    }

    yaxis2_params = {
        "showticklabels": False,
        "showline": True,
        "linewidth": 1,
        "linecolor": apc.All["arcadia:brightgrey"],
        "mirror": True,
    }

    i = 0
    for group in count_features_results["annotation_count"].keys():
        bar = semantic_barchart_plotly(
            count_features_results["annotation_count"], group, colors_dict[group]
        )

        fig.add_trace(bar, row=flattened_indices[i][0], col=flattened_indices[i][1])
        next(fig.select_xaxes(row=flattened_indices[i][0], col=flattened_indices[i][1])).update(
            xaxis_params
        )
        next(fig.select_yaxes(row=flattened_indices[i][0], col=flattened_indices[i][1])).update(
            yaxis_params
        )

        i += 1

        image = wordcloud_image(count_features_results["wordclouds"], group, colors_dict[group])

        fig.add_trace(image, row=flattened_indices[i][0], col=flattened_indices[i][1])
        next(fig.select_xaxes(row=flattened_indices[i][0], col=flattened_indices[i][1])).update(
            xaxis2_params
        )
        next(fig.select_yaxes(row=flattened_indices[i][0], col=flattened_indices[i][1])).update(
            yaxis2_params
        )

        i += 1

    fig.update_layout(
        uniformtext_minsize=10,
        uniformtext_mode="show",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width=550 * n_cols,
        height=300 * n_rows,
        font_family="Arial",
        margin_l=10,
        margin_t=20,
        margin_r=10,
        margin_b=10,
        showlegend=False,
    )

    if savefile is not None:
        fig.write_html(savefile)

    if show:
        fig.show()

    return fig


def main():
    args = parse_args()
    features_file = args.features_file
    agg_col = args.agg_column
    annot_col = args.annot_column
    output_file = args.output
    interactive_file = args.interactive
    exclude_words = args.exclude_words
    analysis_name = args.analysis_name
    colors = apc.Palettes["arcadia:AccentAllOrdered"].colors

    if output_file is not None:
        apc.mpl_setup()

        plot_semantic_analysis(
            features_file=features_file,
            agg_col=agg_col,
            annot_col=annot_col,
            colors=colors,
            savefile=output_file,
            exclude_words=exclude_words,
            analysis_name=analysis_name,
        )

    if interactive_file is not None:
        results = count_features(
            features_file=features_file,
            agg_col=agg_col,
            annot_col=annot_col,
            colors=colors,
        )

        semantic_multiplot_plotly(results, colors, savefile=interactive_file)


if __name__ == "__main__":
    main()
