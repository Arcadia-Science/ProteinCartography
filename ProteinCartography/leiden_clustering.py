#!/usr/bin/env python
import argparse
import scanpy as sc
import pandas as pd
import numpy as np

# only import these functions when using import *
__all__ = ["scanpy_leiden_cluster"]


# parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", required=True, help="Input file path of a similarity matrix."
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output path to a file, usually leiden_features.tsv",
    )
    parser.add_argument(
        "-n",
        "--n-neighbors",
        default="10",
        help="Number of n_neighbors to pass to sc.pp.neighbors().",
    )
    parser.add_argument(
        "-c",
        "--n-pcs",
        default="30",
        help="Number of n_pcs to pass to sc.pp.neighbors().",
    )
    parser.add_argument(
        "-l",
        "--cluster-name",
        default="LeidenCluster",
        help="Name of cluster column. Defaults to 'LeidenCluster'.",
    )
    parser.add_argument(
        "-a",
        "--cluster-abbrev",
        default="LC",
        help="Abbreviation to add as prefix for cluster labels. Defaults to 'LC'.",
    )
    args = parser.parse_args()
    return args


def scanpy_leiden_cluster(
    input_file: str,
    savefile=None,
    n_neighbors=10,
    n_pcs=30,
    cluster_name="LeidenCluster",
    cluster_abbrev="LC",
    **kwargs
):
    """
    Uses Scanpy's Leiden clustering implementation to perform clustering.

    Args:
        input_file (str): path of input distances matrix.
        savefile (str): path of destination file.
        n_neighbors (int): number of neighbors for clustering. Defaults to 10.
        n_pcs (int): number of PCs to use for initial PCA.
        **kwargs are passed to `sc.pp.neighbors()`.
    """
    # Load the data
    adata = sc.read_csv(input_file, delimiter="\t")

    # Run intial PCA
    sc.tl.pca(adata, svd_solver="arpack")

    n_neighbors_recommended = int(np.round(len(adata.var) / 10))
    if n_neighbors_recommended > n_neighbors:
        n_neighbors_used = n_neighbors_recommended
    else:
        n_neighbors_used = n_neighbors

    # Run nearest neighbors, umap, then leiden
    # We should probably determine a good empirical default for this
    sc.pp.neighbors(adata, n_neighbors=n_neighbors_used, n_pcs=n_pcs, **kwargs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    # Extract leiden cluster assignment
    membership = pd.DataFrame(adata.obs["leiden"]).reset_index()
    membership.rename(columns={"index": "protid", "leiden": cluster_name}, inplace=True)
    max_chars = len(str(membership[cluster_name].astype(int).max()))
    membership[cluster_name] = cluster_abbrev + membership[cluster_name].apply(
        lambda x: str(x).zfill(max_chars)
    ).astype(str)

    if savefile is not None:
        membership.to_csv(savefile, sep="\t", index=None)

    return membership


# run this if called from the interpreter
def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    neighbors = int(args.n_neighbors)
    pcs = int(args.n_pcs)
    cluster_name = args.cluster_name
    cluster_abbrev = args.cluster_abbrev

    scanpy_leiden_cluster(
        input_file=input_file,
        savefile=output_file,
        n_neighbors=neighbors,
        n_pcs=pcs,
        cluster_name=cluster_name,
        cluster_abbrev=cluster_abbrev,
    )


# check if called from interpreter
if __name__ == "__main__":
    main()
