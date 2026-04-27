#!/usr/bin/env python3
"""
02_integration_clustering.py
=============================
Merge QC-filtered samples, correct batch effects, and cluster cells.

WHAT THIS SCRIPT DOES:
  1. Load one or more QC-filtered .h5ad files from script 01
  2. Normalize counts and select highly variable genes (per sample)
  3. Merge samples and run PCA
  4. Apply Harmony batch correction
  5. Build k-nearest-neighbour graph and embed in UMAP
  6. Leiden clustering at multiple resolutions
  7. Save integrated AnnData

HOW TO USE:
  1. Edit CONFIGURATION below — especially INPUT_PATHS and BATCH_KEY.
  2. Run:  python 02_integration_clustering.py

BEST PRACTICE:
  Normalize and find HVGs independently per sample before merging.
  Never normalize the merged object globally first.

OUTPUT:
  results/02_integration/adata_integrated.h5ad
  results/02_integration/plots/
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px
import plotly.graph_objects as go

sys.path.insert(0, str(Path("..") / "src"))
from sc_workflow.utils import setup_scanpy, save_adata, pretty_counts, leiden_sweep

warnings.filterwarnings("ignore")
setup_scanpy()


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

# List of .h5ad files produced by script 01 (one per sample).
# Add as many paths as you have samples.
INPUT_PATHS = [
    "results/01_qc/adata_qc.h5ad",
    # "results/01_qc/sample2_qc.h5ad",
]

# Column in adata.obs that identifies which sample each cell came from.
BATCH_KEY = "sample"

# --- Normalisation & HVG parameters ------------------------------------------
# Normalise each cell to this total count before log-transformation.
# "1e4" (10,000) is the standard; "median" uses the median count across cells.
TARGET_SUM = "1e4"   # options: "1e4", "1e5", "median"

# Number of highly variable genes to select per sample.
N_HVGS = 3000

# Number of principal components to compute.
N_PCS = 50

# --- Harmony batch correction -------------------------------------------------
RUN_HARMONY = True
HARMONY_THETA = 2.0     # diversity penalty (higher = stronger correction)
HARMONY_MAX_ITER = 10   # number of Harmony iterations

# --- Neighbourhood graph & clustering ----------------------------------------
N_NEIGHBOURS = 30
UMAP_MIN_DIST = 0.3

# Leiden clustering is run at each of these resolutions.
# Higher resolution = more, smaller clusters.
CLUSTERING_RESOLUTIONS = [0.3, 0.5, 0.8, 1.0]

# The resolution used for display in plots (must be in CLUSTERING_RESOLUTIONS).
DEFAULT_RESOLUTION = 0.5

# --- Output ------------------------------------------------------------------
OUTPUT_PATH = "results/02_integration/adata_integrated.h5ad"


# ==============================================================================
#  STEP 1 — LOAD & MERGE SAMPLES
# ==============================================================================

def load_and_merge():
    adatas = []
    for path_str in INPUT_PATHS:
        p = Path(path_str)
        if not p.exists():
            print(f"File not found: {p}")
            sys.exit(1)
        a = sc.read_h5ad(p)
        if BATCH_KEY not in a.obs.columns:
            a.obs[BATCH_KEY] = p.stem
        adatas.append(a)

    if len(adatas) == 1:
        return adatas[0].copy()
    adata_merged = sc.concat(adatas, join="outer", index_unique="-")
    adata_merged.obs_names_make_unique()
    return adata_merged


# ==============================================================================
#  STEP 2 — NORMALISATION, HVG SELECTION & PCA
# ==============================================================================

def normalise_and_pca(adata_merged):
    adata_pp = adata_merged.copy()

    # Always keep raw counts for DE analysis in script 04
    adata_pp.layers["counts"] = adata_pp.X.copy()

    target_sum = None if TARGET_SUM == "median" else float(TARGET_SUM)
    samples = adata_pp.obs[BATCH_KEY].unique().tolist()

    if len(samples) > 1:
        # Normalise and find HVGs per sample, then take the union of HVGs
        # that appear in at least 2 samples. This avoids bias from any
        # single sample dominating gene selection.
        hvg_dfs = []
        for s in samples:
            a = adata_pp[adata_pp.obs[BATCH_KEY] == s].copy()
            sc.pp.normalize_total(a, target_sum=target_sum)
            sc.pp.log1p(a)
            sc.pp.highly_variable_genes(a, n_top_genes=N_HVGS, flavor="seurat_v3",
                                        layer="counts")
            hvg_dfs.append(a.var[["highly_variable"]])

        hvg_matrix = pd.concat(
            [df["highly_variable"].rename(s) for df, s in zip(hvg_dfs, samples)],
            axis=1,
        )
        shared_hvgs = hvg_matrix.sum(axis=1) >= min(2, len(samples))
        adata_pp.var["highly_variable"] = shared_hvgs

        sc.pp.normalize_total(adata_pp, target_sum=target_sum)
        sc.pp.log1p(adata_pp)
    else:
        sc.pp.normalize_total(adata_pp, target_sum=target_sum)
        sc.pp.log1p(adata_pp)
        sc.pp.highly_variable_genes(adata_pp, n_top_genes=N_HVGS, flavor="seurat_v3",
                                    layer="counts")

    sc.pp.scale(adata_pp, max_value=10)
    sc.tl.pca(adata_pp, n_comps=N_PCS, svd_solver="arpack", random_state=42)
    return adata_pp


def save_elbow_plot(adata_pp, plot_dir):
    var_ratio = adata_pp.uns["pca"]["variance_ratio"]
    cumvar = var_ratio.cumsum() * 100
    fig = px.line(
        x=list(range(1, len(var_ratio) + 1)),
        y=var_ratio * 100,
        labels={"x": "PC", "y": "Variance Explained (%)"},
        title="PCA Elbow Plot",
        template="plotly_white",
    )
    fig.add_scatter(x=list(range(1, len(cumvar) + 1)), y=cumvar,
                    mode="lines", name="Cumulative",
                    line=dict(dash="dash", color="orange"))
    fig.update_layout(height=380)
    fig.write_html(Path(plot_dir) / "01_pca_elbow.html")


# ==============================================================================
#  STEP 3 — HARMONY BATCH CORRECTION
# ==============================================================================

def run_harmony(adata_pp):
    adata_integrated = adata_pp.copy()

    if not RUN_HARMONY:
        return adata_integrated, "X_pca"

    try:
        import harmonypy as hm
        ho = hm.run_harmony(
            adata_integrated.obsm["X_pca"],
            adata_integrated.obs,
            vars_use=BATCH_KEY,
            theta=HARMONY_THETA,
            max_iter_harmony=HARMONY_MAX_ITER,
            random_state=42,
            verbose=False,
        )
        adata_integrated.obsm["X_pca_harmony"] = ho.Z_corr.T
        return adata_integrated, "X_pca_harmony"
    except ImportError:
        print("harmonypy not installed (pip install harmonypy). Using uncorrected PCA.")
        return adata_integrated, "X_pca"


# ==============================================================================
#  STEP 4 — UMAP & LEIDEN CLUSTERING
# ==============================================================================

def embed_and_cluster(adata_integrated, use_rep):
    sc.pp.neighbors(adata_integrated, n_neighbors=N_NEIGHBOURS, n_pcs=N_PCS,
                    use_rep=use_rep, random_state=42)
    sc.tl.umap(adata_integrated, min_dist=UMAP_MIN_DIST, random_state=42)
    leiden_sweep(adata_integrated, resolutions=CLUSTERING_RESOLUTIONS, random_state=42)
    return adata_integrated


def save_umap_plots(adata_integrated, plot_dir):
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    res_key = f"leiden_{DEFAULT_RESOLUTION:.2f}"
    if res_key not in adata_integrated.obs.columns:
        res_key = next(
            (c for c in adata_integrated.obs.columns if c.startswith("leiden_")), None
        )
    if res_key is None:
        return

    umap = adata_integrated.obsm["X_umap"]
    obs = adata_integrated.obs.copy()
    obs["UMAP1"] = umap[:, 0]
    obs["UMAP2"] = umap[:, 1]

    fig_cluster = px.scatter(
        obs, x="UMAP1", y="UMAP2", color=res_key,
        title=f"UMAP — Leiden clusters (res {DEFAULT_RESOLUTION})",
        template="plotly_white", opacity=0.7,
    )
    fig_cluster.update_traces(marker_size=3)
    fig_cluster.update_layout(height=500)
    fig_cluster.write_html(plot_dir / "02_umap_clusters.html")

    fig_batch = px.scatter(
        obs, x="UMAP1", y="UMAP2", color=BATCH_KEY,
        title=f"UMAP — Batch ({BATCH_KEY})",
        template="plotly_white", opacity=0.7,
    )
    fig_batch.update_traces(marker_size=3)
    fig_batch.update_layout(height=500)
    fig_batch.write_html(plot_dir / "03_umap_batch.html")


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    plot_dir = Path(OUTPUT_PATH).parent / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    adata_merged = load_and_merge()
    adata_pp = normalise_and_pca(adata_merged)
    save_elbow_plot(adata_pp, plot_dir)
    adata_integrated, use_rep = run_harmony(adata_pp)
    adata_integrated = embed_and_cluster(adata_integrated, use_rep)
    save_umap_plots(adata_integrated, plot_dir)

    save_adata(adata_integrated, OUTPUT_PATH)
    print(f"Saved → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
