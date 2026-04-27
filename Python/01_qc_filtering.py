#!/usr/bin/env python3
"""
01_qc_filtering.py
==================
Quality control and cell filtering for a single 10x scRNA-seq sample.

WHAT THIS SCRIPT DOES:
  1. Load Cell Ranger output (MTX directory, .h5, or .h5ad file)
  2. Calculate per-cell QC metrics (UMI counts, gene counts, % mitochondrial)
  3. Detect and remove doublets with Scrublet
  4. Filter cells and genes using the thresholds set below
  5. Save a clean AnnData object to results/01_qc/

HOW TO USE:
  1. Edit the CONFIGURATION section below.
  2. Run:  python 01_qc_filtering.py
  3. Run this script once per sample, then merge in script 02.

OUTPUT:
  results/01_qc/adata_qc.h5ad
  results/01_qc/plots/  (QC figures)
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots as psub

sys.path.insert(0, str(Path("..") / "src"))
from sc_workflow.utils import (
    setup_scanpy, load_data, add_qc_metrics,
    mad_outlier, adaptive_qc_filter, save_adata, pretty_counts,
)

warnings.filterwarnings("ignore")
setup_scanpy()


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

# Path to Cell Ranger output directory, a .h5 file, or an existing .h5ad file.
# This is the output from script 00 (data/raw/<sample_id>/outs/filtered_feature_bc_matrix)
DATA_PATH = "data/raw/sample1"

# A short name for this sample. Used in output file names and plot titles.
SAMPLE_NAME = "sample1"

# "human" or "mouse" — controls which mitochondrial gene prefixes are used
SPECIES = "human"

# --- Cell filtering thresholds -----------------------------------------------
# Cells outside these ranges are removed.
MIN_GENES = 500         # minimum number of genes detected per cell
MAX_GENES = 8000        # maximum number of genes (very high = likely doublet)
MIN_COUNTS = 1000       # minimum total UMI counts per cell
MAX_COUNTS = 50000      # maximum total UMI counts per cell
MAX_MT_PCT = 20         # maximum mitochondrial gene % (high = dying cell)

# Genes expressed in fewer than this many cells are removed.
MIN_CELLS_PER_GENE = 3

# --- Doublet detection -------------------------------------------------------
RUN_SCRUBLET = True
EXPECTED_DOUBLET_RATE = 0.06   # 0.8% per 1000 cells loaded (10x rule of thumb)
DOUBLET_THRESHOLD = 0.25       # score above this = predicted doublet

# --- Output ------------------------------------------------------------------
OUTPUT_PATH = "results/01_qc/adata_qc.h5ad"


# ==============================================================================
#  STEP 1 — LOAD DATA
# ==============================================================================

def load_raw_data():
    p = Path(DATA_PATH)
    if not p.exists():
        print(f"Path not found: {p.resolve()}\nSet DATA_PATH in the configuration section.")
        sys.exit(1)
    return load_data(p, sample_name=SAMPLE_NAME)


# ==============================================================================
#  STEP 2 — QC METRICS
# ==============================================================================

def compute_qc_metrics(adata):
    adata = adata.copy()
    add_qc_metrics(adata, species=SPECIES)
    adata.obs["log10_counts"] = np.log10(adata.obs["total_counts"] + 1)
    adata.obs["log10_genes"] = np.log10(adata.obs["n_genes_by_counts"] + 1)

    obs = adata.obs
    print(f"Cells: {adata.n_obs:,}  Genes: {adata.n_vars:,}")
    print(f"Median UMI/cell: {obs['total_counts'].median():.0f}  "
          f"Median genes/cell: {obs['n_genes_by_counts'].median():.0f}  "
          f"Median MT%: {obs['pct_counts_mt'].median():.2f}%")
    return adata


def save_qc_plots(adata, plot_dir):
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)
    obs = adata.obs

    fig = psub.make_subplots(
        rows=1, cols=3,
        subplot_titles=["UMI Counts", "Genes Detected", "Mitochondrial %"],
    )
    colors = ["#4C72B0", "#55A868", "#C44E52"]
    data = [obs["total_counts"], obs["n_genes_by_counts"], obs["pct_counts_mt"]]
    names = ["UMI", "Genes", "MT%"]
    for i, (d, n, c) in enumerate(zip(data, names, colors), start=1):
        fig.add_trace(
            go.Violin(y=d, name=n, box_visible=True, meanline_visible=True,
                      fillcolor=c, line_color="black", opacity=0.75, showlegend=False),
            row=1, col=i,
        )
    fig.update_layout(height=420, template="plotly_white",
                      title_text=f"QC Distributions — {SAMPLE_NAME}")
    fig.write_html(plot_dir / "01_qc_violins.html")

    fig2 = px.scatter(
        obs, x="total_counts", y="n_genes_by_counts",
        color="pct_counts_mt", color_continuous_scale="RdYlBu_r", opacity=0.4,
        labels={"total_counts": "UMI Counts", "n_genes_by_counts": "Genes",
                "pct_counts_mt": "MT %"},
        title=f"UMI vs Genes (colour = MT%) — {SAMPLE_NAME}",
        template="plotly_white",
    )
    fig2.update_traces(marker_size=3)
    fig2.write_html(plot_dir / "02_scatter_qc.html")


# ==============================================================================
#  STEP 3 — DOUBLET DETECTION
# ==============================================================================

def detect_doublets(adata):
    if not RUN_SCRUBLET:
        adata.obs["doublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False
        return adata

    try:
        import scrublet as scr
    except ImportError:
        print("scrublet not installed (pip install scrublet). Skipping doublet removal.")
        adata.obs["doublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False
        return adata

    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    scrub = scr.Scrublet(X, expected_doublet_rate=EXPECTED_DOUBLET_RATE, random_state=42)
    scores, _ = scrub.scrub_doublets(min_counts=2, min_cells=3, n_prin_comps=30)
    predicted = scores >= DOUBLET_THRESHOLD

    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = predicted

    n = int(predicted.sum())
    print(f"Doublets detected: {n:,} ({100*n/len(predicted):.1f}% at threshold {DOUBLET_THRESHOLD})")
    return adata


# ==============================================================================
#  STEP 4 — FILTER CELLS & GENES
# ==============================================================================

def filter_cells_and_genes(adata):
    n_before = adata.n_obs
    obs = adata.obs

    mask = (
        (obs["n_genes_by_counts"] >= MIN_GENES)
        & (obs["n_genes_by_counts"] <= MAX_GENES)
        & (obs["total_counts"] >= MIN_COUNTS)
        & (obs["total_counts"] <= MAX_COUNTS)
        & (obs["pct_counts_mt"] <= MAX_MT_PCT)
    )
    if RUN_SCRUBLET and "predicted_doublet" in obs.columns:
        mask &= ~obs["predicted_doublet"]

    adata_filtered = adata[mask].copy()

    g_before = adata_filtered.n_vars
    sc.pp.filter_genes(adata_filtered, min_cells=MIN_CELLS_PER_GENE)

    print(f"Cells: {n_before:,} → {adata_filtered.n_obs:,} "
          f"({n_before - adata_filtered.n_obs:,} removed)")
    print(f"Genes: {g_before:,} → {adata_filtered.n_vars:,} "
          f"(removed genes in < {MIN_CELLS_PER_GENE} cells)")
    return adata_filtered


def save_before_after_plots(adata_raw, adata_filtered, plot_dir):
    plot_dir = Path(plot_dir)
    metrics = [
        ("total_counts", "UMI Counts"),
        ("n_genes_by_counts", "Genes Detected"),
        ("pct_counts_mt", "MT %"),
    ]
    fig = psub.make_subplots(rows=1, cols=3,
                              subplot_titles=[m[1] for m in metrics])
    for j, (col, label) in enumerate(metrics, start=1):
        before = adata_raw.obs[col] if col in adata_raw.obs.columns else []
        after = adata_filtered.obs[col]
        fig.add_trace(
            go.Violin(y=before, name="Before", fillcolor="#AEC6E8",
                      line_color="#4C72B0", opacity=0.7, showlegend=(j == 1)),
            row=1, col=j,
        )
        fig.add_trace(
            go.Violin(y=after, name="After", fillcolor="#A8D5A2",
                      line_color="#228B22", opacity=0.7, showlegend=(j == 1)),
            row=1, col=j,
        )
    fig.update_layout(height=420, template="plotly_white",
                      title_text="QC Metrics: Before vs After Filtering",
                      violinmode="overlay")
    fig.write_html(plot_dir / "03_before_after_filter.html")


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    plot_dir = Path(OUTPUT_PATH).parent / "plots"

    adata_raw = load_raw_data()
    adata_metrics = compute_qc_metrics(adata_raw)
    save_qc_plots(adata_metrics, plot_dir)
    adata_scrublet = detect_doublets(adata_metrics)
    adata_filtered = filter_cells_and_genes(adata_scrublet)
    save_before_after_plots(adata_metrics, adata_filtered, plot_dir)

    save_adata(adata_filtered, OUTPUT_PATH)
    print(f"Saved → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
