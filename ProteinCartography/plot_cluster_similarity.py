#!/usr/bin/env python
import argparse

import numpy as np
import pandas as pd
import plotly.express as px
from color_utils import arcadia_viridis

__all__ = ["calculate_group_similarity", "plot_group_similarity"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--matrix-file",
        required=True,
        help="Path to input all-v-all similarity matrix.",
    )
    parser.add_argument(
        "-f",
        "--features-file",
        required=True,
        help="Path to features file for grouping.",
    )
    parser.add_argument(
        "-c", "--features-column", required=True, help="Column to aggregate groups on."
    )
    parser.add_argument("-T", "--output-tsv", help="Path to output TSV file.")
    parser.add_argument("-H", "--output-html", help="Path to output HTML file.")
    args = parser.parse_args()

    return args


def calculate_group_similarity(
    matrix_file: str, features_file: str, features_column: str, output_file=None
):
    """
    Takes an all-v-all similarity matrix_file and averages across a grouping
    found in a features_file.
    The grouping is specified based on a column of the features file.

    Each row and column in the matrix_file should be a protid.
    Within the features file, each protid should have an associated value
    in the specified features_column, which is its category.
    The similarity matrix will be aggregated based on those groupings,
    and the similarity values averaged across the whole group versus every other group.

    Args:
        matrix_file (str): path of input similarity matrix file.
        features_file (str): path of input features file.
        features_column (str); column of the features file to aggregate on.
        output_file (str): path of destination file.
    """
    # load in each dataframe
    pivot_df = pd.read_csv(matrix_file, sep="\t")
    features_df = pd.read_csv(features_file, sep="\t")

    # merge the two dataframes on protid
    pivot_agg = pivot_df.merge(features_df, on="protid")

    # group entries along one axis by feature groups, applying mean
    pivot_agg = pivot_agg.groupby(features_column).agg(
        {i: list if i == "protid" or i == features_column else np.mean for i in pivot_agg.columns}
    )

    # tidy up dataframe
    pivot_agg.reset_index(drop=True, inplace=True)
    pivot_agg.drop(columns=["protid", features_column], inplace=True)

    # transpose dataframe
    pivot_t = pivot_agg.transpose()
    pivot_t = pivot_t.reset_index().rename(columns={"index": "protid"})

    # merge features dataframe again and groupby feature
    pivot_t_agg = pivot_t.merge(features_df, on="protid")
    pivot_t_agg = pivot_t_agg.groupby(features_column).agg(
        {i: list if i == "protid" or i == features_column else np.mean for i in pivot_t_agg.columns}
    )

    # clean up and reset groupings
    pivot_t_agg.drop(columns=["protid", features_column], inplace=True)
    pivot_t_agg.columns = pivot_t_agg.index

    if output_file is not None:
        pivot_t_agg.to_csv(output_file, sep="\t")

    return pivot_t_agg


def plot_group_similarity(
    group_similarity_file: str,
    output_file=None,
    plot_width=700,
    plot_height=700,
    show=False,
):
    """
    Takes an all-v-all aggregated group_similarity_file and makes a plotly interactive heatmap.
    Saves to an HTML file if an output filepath is provided.

    Args:
        group_similarity_file (str): path of input group-averaged similarity matrix file.
        output_file (str): path of destination file.
        plot_width (int): width of plot in pixels. Defaults to 700.
        plot_height (int): height of plot in pixels. Defaults to 700.
        show (bool): whether or not to show the plot.

    Returns:
        A Plotly figure object for the corresponding interactive plot.
    """
    # opens the aggregated similarity matrix file
    sim_df = pd.read_csv(group_similarity_file, index_col=0, sep="\t")

    # generates a plotly heatmap for that matrix file
    fig = px.imshow(sim_df, color_continuous_scale=arcadia_viridis, range_color=[0, 1])

    colorbar_dict = dict(
        title="Similarity",
        x=1,
        y=0.5,
        xanchor="left",
        yanchor="middle",
        title_font_size=14,
        title_side="top",
        len=0.7,
        thickness=20,
    )

    fig.update_layout(width=plot_width, height=plot_height, coloraxis_colorbar=colorbar_dict)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0), paper_bgcolor="rgba(0,0,0,0)")

    fig.update_xaxes(side="top", title="Target")
    fig.update_yaxes(title="Query")

    try:
        fig.update_layout(font=dict(family="Arial"))
    except Exception:
        pass

    plot_config = {
        "displayModeBar": True,
        "toImageButtonOptions": {"filename": "heatmap", "format": "svg"},
        "modeBarButtonsToRemove": ["zoomIn", "zoomOut"],
    }

    if output_file is not None:
        fig.write_html(output_file, config=plot_config)

    if show:
        fig.show(config=plot_config)

    return fig


def main():
    args = parse_args()
    matrix = args.matrix_file
    features = args.features_file
    column = args.features_column
    output_tsv = args.output_tsv
    output_html = args.output_html

    calculate_group_similarity(matrix, features, column, output_file=output_tsv)

    if output_tsv is not None:
        plot_group_similarity(output_tsv, output_file=output_html)


if __name__ == "__main__":
    main()
