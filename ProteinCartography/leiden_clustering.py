#!/usr/bin/env python
"""Leiden clustering over a ProteinCartography TM-score matrix.

For N < 3, assign every protein to a single cluster and skip Scanpy. Otherwise
clamp ``n_neighbors`` / ``n_pcs`` to valid ranges for the matrix size.
"""

from __future__ import annotations
import argparse

import numpy as np
import pandas as pd
import scanpy as sc

__all__ = ["scanpy_leiden_cluster"]


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


def _singleton_membership(
    protids: list[str], cluster_name: str, cluster_abbrev: str
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protid": protids,
            cluster_name: [f"{cluster_abbrev}0"] * len(protids),
        }
    )


def scanpy_leiden_cluster(
    input_file: str,
    savefile=None,
    n_neighbors=10,
    n_pcs=30,
    cluster_name="LeidenCluster",
    cluster_abbrev="LC",
    **kwargs,
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
    adata = sc.read_csv(input_file, delimiter="\t")
    n = int(adata.n_obs)
    if n < 3:
        print(
            f"[leiden_clustering] N={n} too small for Leiden; "
            f"assigning all proteins to {cluster_abbrev}0"
        )
        membership = _singleton_membership(list(adata.obs_names), cluster_name, cluster_abbrev)
        if savefile is not None:
            membership.to_csv(savefile, sep="\t", index=None)
        return membership

    sc.tl.pca(adata, svd_solver="arpack")

    n_neighbors = int(n_neighbors)
    n_pcs = int(n_pcs)
    n_neighbors_recommended = int(np.round(n / 10))
    n_neighbors_used = max(n_neighbors, n_neighbors_recommended)
    # Scanpy requires 1 < n_neighbors <= N - 1 for a connected graph.
    n_neighbors_used = max(2, min(n_neighbors_used, n - 1))
    max_pcs = max(1, min(n - 1, adata.n_vars - 1 if adata.n_vars > 1 else 1))
    n_pcs_used = max(1, min(n_pcs, max_pcs))

    sc.pp.neighbors(adata, n_neighbors=n_neighbors_used, n_pcs=n_pcs_used, **kwargs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    membership = pd.DataFrame(adata.obs["leiden"]).reset_index()
    membership.rename(columns={"index": "protid", "leiden": cluster_name}, inplace=True)
    max_chars = len(str(membership[cluster_name].astype(int).max()))
    membership[cluster_name] = cluster_abbrev + membership[cluster_name].apply(
        lambda x: str(x).zfill(max_chars)
    ).astype(str)

    if savefile is not None:
        membership.to_csv(savefile, sep="\t", index=None)

    return membership


def main():
    args = parse_args()
    n_neighbors = int(args.n_neighbors)
    n_pcs = int(args.n_pcs)

    scanpy_leiden_cluster(
        args.input,
        args.output,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        cluster_name=args.cluster_name,
        cluster_abbrev=args.cluster_abbrev,
    )


if __name__ == "__main__":
    main()
