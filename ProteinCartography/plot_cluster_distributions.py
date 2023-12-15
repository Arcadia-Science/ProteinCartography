#!/usr/bin/env python

# import any necessary functions
# argparse is always required
import argparse
from itertools import chain

import arcadia_pycolor as apc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stat

apc.mpl_setup()

# only import these functions when using from template import *
# fill this in with any functions that are not "parse_args()" or "main()"
__all__ = ["remove_nans", "distribution_test", "plot_distribution_violins"]

DEFAULT_PLOTTING_RULES = {
    "Length": {
        "textlabel": "Length",
        "facecolor": apc.All["arcadia:oat"],
        "edgecolor": apc.All["arcadia:canary"],
    },
    "pdb_confidence": {
        "textlabel": "Mean pLDDT",
        "facecolor": apc.All["arcadia:bluesky"],
        "edgecolor": apc.All["arcadia:aegean"],
    },
    "Annotation": {
        "textlabel": "Annotation Score",
        "facecolor": apc.All["arcadia:wish"],
        "edgecolor": apc.All["arcadia:aster"],
    },
}


# parse command line arguments
def parse_args():
    # initialize an argument parser object
    parser = argparse.ArgumentParser()

    # add arguments one at a time to that argument parser object
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to input features file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        nargs="?",
        default="",
        help="Path to output plot file.",
    )
    parser.add_argument(
        "-k",
        "--keyid",
        help="keyid for plotting. Usually input protid.",
    )

    # actually run the argument parsing and accept input from the command line
    args = parser.parse_args()

    # return valid arguments as an object
    return args


def remove_nans(nested_list: list, default=(0,)):
    """
    Remove NaNs from a nested list.
    Returns (0,) if the list is empty after removing NaNs.
    """
    cleaned = [
        remove_nans(item) if isinstance(item, list) else item
        for item in nested_list
        if not (isinstance(item, float) and pd.isna(item))
    ]
    return default if not cleaned else cleaned


def distribution_test(
    data_dict: dict,
    data1_keys: list,
    data2_keys: list = None,
    method: str = "MWU",
    *args,
    **kwargs,
):
    """
    Perform a statistical test to compare the distributions of two sets of data.

    Args:
    data_dict: A dictionary containing the data to be compared.
        The keys of the dictionary are the cluster IDs.
        The values of the dictionary are lists of the distribution of values per cluster.
    data1_keys: The keys of the data to be aggregated into
        the first distribution for comparison.
    data2_keys: The keys of the data to be aggregated into
        the second distribution for comparison.
        If None, all keys not in data1_keys are used.
    method: The statistical test to be used. Options are "MWU" (Mann-Whitney U test)
        and "KS2S" (Kolmogorov-Smirnov 2-sample test).
    *args: Additional arguments to be passed to the statistical test function.
    **kwargs: Additional keyword arguments to be passed to the statistical test function.
    """

    if any([key not in data_dict.keys() for key in data1_keys]):
        raise ValueError("Invalid data1_keys")
    if data2_keys is not None:
        if any([key not in data_dict.keys() for key in data2_keys]):
            raise ValueError("Invalid data2_keys")

    data1 = []

    # aggregate data for data1 into a single flat list
    data1 = list(chain.from_iterable([data_dict[k] for k in data1_keys]))

    data2 = []

    # if data2_keys is None, use all keys not in data1_keys
    if data2_keys is None:
        data2_keys = [key for key in data_dict.keys() if key not in data1_keys]

    # aggregate data for data2 into a single flat list
    data2 = list(chain.from_iterable([data_dict[k] for k in data2_keys]))

    if method == "MWU":
        method_fxn = stat.mannwhitneyu
    elif method == "KS2S":
        method_fxn = stat.ks_2samp
    else:
        raise ValueError(f"Invalid method: {method}")

    result = method_fxn(data1, data2, *args, **kwargs)

    return result


def generate_violin_plotting_rules(keyid: str = None):
    """
    Generate the plotting rules for the violin plots.

    Args:
    keyid (str): The keyid for the input protein to use for identifying TMscore and fident columns.
    """
    plotting_rules = DEFAULT_PLOTTING_RULES

    if keyid is not None:
        plotting_rules[f"TMscore_v_{keyid}"] = {
            "textlabel": f"TMscore v {keyid}",
            "facecolor": apc.All["arcadia:mint"],
            "edgecolor": apc.All["arcadia:seaweed"],
        }
        plotting_rules[f"fident_v_{keyid}"] = {
            "textlabel": f"fident v {keyid}",
            "facecolor": apc.All["arcadia:blossom"],
            "edgecolor": apc.All["arcadia:rose"],
        }
    return plotting_rules


def plot_distribution_violins(
    input_filepath: str,
    plotting_rules: dict,
    cluster_column: str = "LeidenCluster",
    output_filepath: str = None,
    threshold_padj_value: float = 0.001,
    target_clusters: list = None,
    target_cluster_label: str = None,
    keyid: str = None,
):
    """
    Plot violin plots of the distributions of the given columns for each cluster.

    Args:
    input_filepath: The filepath to the input aggregated_features file containing
        at least a "protid" column and a "LeidenCluster" column.
    plotting_rules: A dictionary containing the plotting rules for each column to be plotted.
        For each column, the plotting rules should be a dictionary containing the following keys:
        "textlabel" (str): The label for the x-axis.
        "facecolor" (str): The facecolor for the violin plot.
        "edgecolor" (str): The edgecolor for the violin plot.
    cluster_column: The column name of the cluster column. Defaults to "LeidenCluster".
    output_filepath: The filepath to save the output plot.
    threshold_padj_value: The threshold p-value for the corrected statistical test.
    target_clusters: The list of clusters to be used for comparison to all other clusters.
        If None, the cluster containing the "keyid" protein will be used as the target cluster.
    target_cluster_label: The label for the input group in the plot.
    keyid: The keyid for the input protein to use for identifying a target cluster.
    """
    combined_data_df = pd.read_csv(input_filepath, sep="\t")
    grouped = combined_data_df.groupby(cluster_column).agg(list)

    if keyid is not None and target_clusters is None:
        target_cluster_row = combined_data_df[combined_data_df["protid"] == keyid]
        target_clusters = list(target_cluster_row[cluster_column].unique())
        target_cluster_label = f"input {target_clusters[0]}"

    # Create a dictionary to store the distributions of the data
    distribution_dict = {}

    # ignore missing plotting rule columns
    all_columns = set(plotting_rules.keys())
    valid_columns = all_columns.intersection(grouped.columns)
    invalid_columns = all_columns.difference(grouped.columns)
    if len(invalid_columns):
        print(f"Ignoring the following invalid plotting rule columns: {invalid_columns}.")

    for col in valid_columns:
        distribution_dict[col] = dict(grouped[col])
        for item in distribution_dict[col]:
            distribution_dict[col][item] = remove_nans(distribution_dict[col][item])

    # Create a dictionary to store the uncorrected statistical test results
    stat_results_dict = {}

    for col in distribution_dict:
        stat_results_dict[col] = {
            cluster_id: distribution_test(distribution_dict[col], [cluster_id], target_clusters)
            for cluster_id in distribution_dict[col]
        }

    # Create a dictionary to store the corrected statistical test results
    corrected_stat_results_dict = {}

    for col in stat_results_dict:
        pvals = [stat_results_dict[col][cluster_id].pvalue for cluster_id in stat_results_dict[col]]
        corrected_pvals = stat.false_discovery_control(pvals)
        corrected_stat_results_dict[col] = {
            cluster_id: corrected_pvals[i]
            for i, cluster_id in enumerate(stat_results_dict[col].keys())
        }

    fig, axs = plt.subplots(
        nrows=1,
        ncols=len(valid_columns),
        figsize=(1.4 * len(valid_columns), 6.5),
    )

    # generate plots by iterating through plotting rules columns
    for i, (col, ax) in enumerate(zip(valid_columns, axs)):
        values = list(distribution_dict[col].values())

        if target_clusters is not None:
            # Aggregate data for the target clusters
            target_clusters_dist = []
            for target_cluster in target_clusters:
                target_clusters_dist.extend(distribution_dict[col][target_cluster])
            # Calculate median of the target clusters
            target_clusters_median = np.median(target_clusters_dist)

            values.append(target_clusters_dist)

            ytick_increment = 2
            ylim_increment = 1.75
        else:
            ytick_increment = 1
            ylim_increment = 0.75

        # calculate medians for each cluster
        medians_dict = {
            cluster_id: np.median(values) for cluster_id, values in distribution_dict[col].items()
        }
        parts = ax.violinplot(values, vert=False, showextrema=False, widths=0.8)

        for pc in parts["bodies"]:
            pc.set_facecolor(plotting_rules[col]["facecolor"])
            pc.set_edgecolor(plotting_rules[col]["edgecolor"])
            pc.set_linewidths(1)
            pc.set_alpha(0.8)

        for j, (cluster_id, median) in enumerate(medians_dict.items()):
            if corrected_stat_results_dict[col][cluster_id] > threshold_padj_value:
                marker = "."
                color = plotting_rules[col]["edgecolor"]
            else:
                marker = "*"
                color = plotting_rules[col]["edgecolor"]

            if target_clusters is not None:
                if corrected_stat_results_dict[col][cluster_id] > threshold_padj_value:
                    marker = "."
                    color = plotting_rules[col]["edgecolor"]

                elif median < target_clusters_median:
                    marker = "v"
                    color = apc.All["arcadia:marineblue"]

                elif median > target_clusters_median:
                    marker = "^"
                    color = apc.All["arcadia:crow"]

                else:
                    marker = "o"
                    color = plotting_rules[col]["edgecolor"]

            ax.scatter(median, j + 1, marker=marker, color=color, s=80)

        # add an additional marker for the target cluster distribution
        if target_clusters is not None:
            marker = "."
            color = plotting_rules[col]["edgecolor"]
            ax.scatter(
                target_clusters_median, len(medians_dict) + 1, marker=marker, color=color, s=80
            )

        ax.set_xlabel(plotting_rules[col]["textlabel"])

        # only plot y-axis labels on the first plot
        if i == 0:
            if target_cluster_label is not None:
                tick_labels = list(distribution_dict[col].keys()) + [target_cluster_label]
            else:
                tick_labels = list(distribution_dict[col].keys())
            ax.set_yticks(
                range(1, len(distribution_dict[col].keys()) + ytick_increment),
                tick_labels,
            )
        else:
            ax.set_yticks([])
        ax.set_ylim((0.25, len(distribution_dict[col].keys()) + ylim_increment))

    plt.subplots_adjust(wspace=0.1)
    plt.tight_layout()

    if output_filepath is not None:
        plt.savefig(output_filepath, bbox_inches="tight")


# where the actual functions you've defined get called
# this part of the script is actually executed
def main():
    # run arg parsing function and return arguments as object
    args = parse_args()

    # access each of the arguments as a variable
    input_file = args.input
    output_file = args.output
    keyid = args.keyid

    plotting_rules = generate_violin_plotting_rules(keyid)

    plot_distribution_violins(
        input_filepath=input_file,
        output_filepath=output_file,
        plotting_rules=plotting_rules,
        keyid=keyid,
    )


# this checks whether the script is called from the command line
# if it is, then runs the main() function
if __name__ == "__main__":
    main()
