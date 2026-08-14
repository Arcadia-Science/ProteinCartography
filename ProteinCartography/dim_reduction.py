#!/usr/bin/env python
"""Dimensionality reduction for ProteinCartography similarity matrices.

Clamps UMAP/t-SNE parameters for small N and falls back to a PCA/identity
layout when N < 3 so tiny maps still produce coordinates.
"""

from __future__ import annotations
import argparse

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

__all__ = ["calculate_PCA", "calculate_TSNE", "calculate_UMAP"]

MODES = ["pca", "tsne", "umap", "pca_tsne", "pca_umap"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-p", "--output-prefix", help="Prefix for resulting .tsv files.")
    parser.add_argument(
        "-m",
        "--mode",
        default="pca",
        help=f"Mode of dimensionality reduction.\nValid arguments are {MODES}.",
    )
    parser.add_argument(
        "-r",
        "--random-state",
        default="123456",
        help="Random state for umap and tsne modes.",
    )
    return parser.parse_args()


def _save_path(pivot_file: str, saveprefix: str | None, dimtype: str) -> str:
    if saveprefix is not None:
        return "_".join([saveprefix, dimtype + ".tsv"])
    return pivot_file.replace(".tsv", "_" + dimtype + ".tsv")


def calculate_PCA(
    pivot_file: str,
    random_state: int,
    n_components=2,
    save=False,
    saveprefix=None,
    dimtype="pca",
    prep_step=False,
    **kwargs,
):
    pivoted_df = pd.read_csv(pivot_file, sep="\t", index_col="protid")

    max_n_components = min(pivoted_df.shape)
    if max_n_components < 1:
        raise ValueError(f"empty similarity matrix: {pivot_file}")
    if n_components > max_n_components:
        print(
            f"Warning: the specified value of `n_components` ({n_components})"
            f"cannot be greater than the number of samples ({max_n_components}),"
            f"so `n_components` will be set to {max_n_components}."
        )
        n_components = max_n_components

    # `svd_solver` is set explicitly because the default ('auto') switches to the randomized
    # solver once the matrix has more than 500 rows or columns, and the randomized solver draws
    # its projection from the global numpy random state. Every run of a map with more than ~500
    # proteins therefore produced different principal components, which are then the input to
    # t-SNE and UMAP, so seeding those alone did not make the pipeline reproducible.
    # 'full' is exact and deterministic, and is also invariant to the order of the columns.
    pca = PCA(n_components=n_components, svd_solver="full", random_state=random_state, **kwargs)
    pca_results = pca.fit_transform(pivoted_df)
    pca_results_df = pd.DataFrame(
        pca_results,
        columns=[f"PC{i}" for i in range(pca_results.shape[1])],
        index=pivoted_df.index,
    )

    savefile = _save_path(pivot_file, saveprefix, dimtype)
    if save:
        pca_results_df.to_csv(savefile, sep="\t")
    if prep_step:
        return savefile
    return pca_results_df


def calculate_TSNE(
    pivot_file: str,
    random_state: int,
    n_components=2,
    perplexity=50,
    n_iter=2000,
    save=False,
    saveprefix=None,
    dimtype="tsne",
    **kwargs,
):
    pivoted_df = pd.read_csv(pivot_file, sep="\t", index_col="protid")
    n = len(pivoted_df)
    if n < 2:
        # Single point — emit a degenerate layout.
        coords = np.zeros((n, n_components), dtype=float)
        tsne_results_df = pd.DataFrame(
            coords,
            columns=[f"tSNE{i + 1}" for i in range(n_components)],
            index=pivoted_df.index,
        )
    else:
        # sklearn requires perplexity < n_samples. Clamp only for that;
        # do not also scale perplexity for every N under 250 (that shifts
        # ordinary maps, e.g. N=100 from 50 → 21).
        perplexity_check = min(perplexity, max(1, n - 1))
        perplexity_check = min(perplexity_check, n - 1)
        tsne = TSNE(
            n_components=min(n_components, max(1, n - 1)),
            perplexity=float(perplexity_check),
            n_iter=n_iter,
            random_state=random_state,
            **kwargs,
        )
        tsne_results = tsne.fit_transform(pivoted_df)
        # Pad to requested components if sklearn reduced them.
        if tsne_results.shape[1] < n_components:
            pad = np.zeros((n, n_components - tsne_results.shape[1]))
            tsne_results = np.hstack([tsne_results, pad])
        tsne_results_df = pd.DataFrame(
            tsne_results[:, :n_components],
            columns=[f"tSNE{i + 1}" for i in range(n_components)],
            index=pivoted_df.index,
        )

    savefile = _save_path(pivot_file, saveprefix, dimtype)
    if save:
        tsne_results_df.to_csv(savefile, sep="\t")
    return tsne_results_df


def _small_n_layout(pivoted_df: pd.DataFrame, n_components: int, prefix: str) -> pd.DataFrame:
    """PCA or identity layout when N is too small for UMAP/t-SNE."""
    n = len(pivoted_df)
    if n == 0:
        raise ValueError("empty matrix")
    if n == 1:
        coords = np.zeros((1, n_components), dtype=float)
    else:
        # Prefer existing PC columns if this is already a PCA matrix.
        pc_cols = [c for c in pivoted_df.columns if str(c).startswith("PC")]
        if len(pc_cols) >= n_components:
            coords = pivoted_df[pc_cols[:n_components]].to_numpy(dtype=float)
        else:
            k = min(n_components, n, pivoted_df.shape[1])
            pca = PCA(n_components=k)
            reduced = pca.fit_transform(pivoted_df)
            if reduced.shape[1] < n_components:
                pad = np.zeros((n, n_components - reduced.shape[1]))
                coords = np.hstack([reduced, pad])
            else:
                coords = reduced[:, :n_components]
    return pd.DataFrame(
        coords,
        columns=[f"{prefix}{i + 1}" for i in range(n_components)],
        index=pivoted_df.index,
    )


def calculate_UMAP(
    pivot_file: str,
    random_state: int,
    n_components=2,
    n_neighbors=80,
    min_dist=0.5,
    save=False,
    saveprefix=None,
    dimtype="umap",
    **kwargs,
):
    pivoted_df = pd.read_csv(pivot_file, sep="\t", index_col="protid")
    n = len(pivoted_df)

    # umap-learn requires n_neighbors >= 2 and n_neighbors < n (effectively ≤ n-1).
    if n < 3:
        print(
            f"[dim_reduction] N={n} too small for UMAP; writing PCA/identity layout "
            f"as {dimtype} coordinates",
            flush=True,
        )
        umap_results_df = _small_n_layout(pivoted_df, n_components, "UMAP")
    else:
        neighbors_check = min(int(n_neighbors), n - 1)
        neighbors_check = max(2, neighbors_check)
        umap_fxn = UMAP(
            n_components=n_components,
            random_state=random_state,
            n_neighbors=neighbors_check,
            min_dist=min_dist,
            **kwargs,
        )
        umap_results = umap_fxn.fit_transform(pivoted_df)
        umap_results_df = pd.DataFrame(
            umap_results,
            columns=[f"UMAP{i + 1}" for i in range(umap_results.shape[1])],
            index=pivoted_df.index,
        )

    savefile = _save_path(pivot_file, saveprefix, dimtype)
    if save:
        umap_results_df.to_csv(savefile, sep="\t")
    return umap_results_df


def main():
    args = parse_args()
    pivot_file = args.input
    saveprefix = args.output_prefix
    mode = args.mode.lower()

    try:
        random_state = int(args.random_state)
    except Exception:
        random_state = 123456

    if mode not in MODES:
        raise Exception(f"{mode} provided is not valid.\nValid modes include {MODES}.")

    if mode == "pca":
        calculate_PCA(pivot_file, random_state, save=True, saveprefix=saveprefix)
    elif mode == "tsne":
        calculate_TSNE(pivot_file, random_state, save=True, saveprefix=saveprefix)
    elif mode == "umap":
        calculate_UMAP(pivot_file, random_state, save=True, saveprefix=saveprefix)
    elif mode == "pca_tsne":
        saveprefix1 = pivot_file.replace(".tsv", "temp1")
        pca_results_file = calculate_PCA(
            pivot_file,
            random_state,
            save=True,
            saveprefix=saveprefix1,
            n_components=30,
            prep_step=True,
        )
        saveprefix2 = pca_results_file.replace("temp1", "").replace(".tsv", "")
        calculate_TSNE(pca_results_file, random_state, save=True, saveprefix=saveprefix2)
    elif mode == "pca_umap":
        saveprefix1 = pivot_file.replace(".tsv", "temp2")
        pca_results_file = calculate_PCA(
            pivot_file,
            random_state,
            save=True,
            saveprefix=saveprefix1,
            n_components=30,
            prep_step=True,
        )
        saveprefix2 = pca_results_file.replace("temp2", "").replace(".tsv", "")
        calculate_UMAP(pca_results_file, random_state, save=True, saveprefix=saveprefix2)


if __name__ == "__main__":
    main()
