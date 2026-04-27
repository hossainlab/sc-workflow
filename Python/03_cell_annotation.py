#!/usr/bin/env python3
"""
03_cell_annotation.py
=====================
Assign cell type identities to Leiden clusters.

WHAT THIS SCRIPT DOES:
  1. CellTypist — automated reference-based annotation (pre-trained models)
  2. Marker gene scoring — score each cell against canonical gene sets
  3. Cluster marker discovery — Wilcoxon test to find top DE genes per cluster
  4. Print a per-cluster consensus table so you can review and finalise labels
  5. Apply the final annotation to adata.obs["cell_type"]
  6. Save the annotated AnnData

HOW TO USE:
  1. Edit CONFIGURATION below.
  2. Run:  python 03_cell_annotation.py
  3. Review the CLUSTER ANNOTATION TABLE printed to the terminal.
  4. Edit the MANUAL_ANNOTATION dict to override any automated label.
  5. Re-run to apply your manual labels and save.

OUTPUT:
  results/03_annotation/adata_annotated.h5ad
  results/03_annotation/plots/
  results/03_annotation/cluster_markers.csv
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px

sys.path.insert(0, str(Path("..") / "src"))
from sc_workflow.utils import setup_scanpy, save_adata, pretty_counts, get_marker_genes

warnings.filterwarnings("ignore")
setup_scanpy()


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH = "results/02_integration/adata_integrated.h5ad"

# The Leiden cluster column to annotate (produced by script 02).
CLUSTER_KEY = "leiden_0.50"

# "human" or "mouse"
ORGANISM = "human"

# --- CellTypist ---------------------------------------------------------------
# Choose a model appropriate for your tissue type.
# Available models: https://www.celltypist.org/models
CELLTYPIST_MODEL = "Immune_All_High.pkl"

# If True, smooths predictions within clusters (recommended).
CELLTYPIST_MAJORITY_VOTING = True

# --- Cluster markers (Wilcoxon test) -----------------------------------------
# Number of top marker genes to find per cluster.
N_TOP_MARKER_GENES = 10

# --- Manual annotation override -----------------------------------------------
# After reviewing the printed table, add entries here to override automated labels.
# Format:  "cluster_id": "Your Cell Type Name"
# Example: {"0": "CD4+ T cell", "3": "Monocyte"}
# Leave empty ({}) to use automated CellTypist labels for all clusters.
MANUAL_ANNOTATION = {}

# --- Output ------------------------------------------------------------------
OUTPUT_PATH = "results/03_annotation/adata_annotated.h5ad"


# ==============================================================================
#  STEP 1 — LOAD DATA
# ==============================================================================

def load_data():
    p = Path(INPUT_PATH)
    if not p.exists():
        print(f"File not found: {p}")
        sys.exit(1)
    adata = sc.read_h5ad(p)

    if CLUSTER_KEY not in adata.obs.columns:
        available = [c for c in adata.obs.columns if c.startswith("leiden_")]
        print(f"Column '{CLUSTER_KEY}' not found. Available: {available}")
        sys.exit(1)

    return adata


# ==============================================================================
#  STEP 2 — CELLTYPIST ANNOTATION
# ==============================================================================

def run_celltypist(adata):
    try:
        import celltypist
        from celltypist import models

        adata_ct = adata.copy()
        sc.pp.normalize_total(adata_ct, target_sum=1e4)
        sc.pp.log1p(adata_ct)

        model = models.Model.load(model=CELLTYPIST_MODEL)
        preds = celltypist.annotate(
            adata_ct,
            model=model,
            majority_voting=CELLTYPIST_MAJORITY_VOTING,
            over_clustering=CLUSTER_KEY,
        )
        adata_out = preds.to_adata()

        label_col = "majority_voting" if CELLTYPIST_MAJORITY_VOTING else "predicted_labels"
        adata.obs["celltypist_label"] = adata_out.obs[label_col]
        adata.obs["celltypist_conf"] = adata_out.obs["conf_score"].values

    except ImportError:
        print("celltypist not installed (pip install celltypist). Using 'unknown' labels.")
        adata.obs["celltypist_label"] = "unknown"
        adata.obs["celltypist_conf"] = 0.0

    return adata


# ==============================================================================
#  STEP 3 — MARKER GENE SCORING
# ==============================================================================

def run_marker_scoring(adata):
    markers = get_marker_genes(ORGANISM)
    scores = {}

    for ct, genes in markers.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            continue
        key = f"score_{ct.replace(' ', '_').replace('+', 'pos').replace('/', '_')}"
        sc.tl.score_genes(adata, gene_list=present, score_name=key)
        scores[ct] = key

    if scores:
        score_matrix = adata.obs[list(scores.values())].values
        best_idx = np.argmax(score_matrix, axis=1)
        ct_names = list(scores.keys())
        adata.obs["marker_score_label"] = [ct_names[i] for i in best_idx]
    else:
        adata.obs["marker_score_label"] = "unknown"
        print("Warning: no marker genes overlapped with dataset. Check ORGANISM setting.")

    return adata


# ==============================================================================
#  STEP 4 — DATA-DRIVEN CLUSTER MARKERS
# ==============================================================================

def find_cluster_markers(adata, out_dir):
    sc.tl.rank_genes_groups(
        adata,
        groupby=CLUSTER_KEY,
        method="wilcoxon",
        use_raw=False,
        key_added="rank_genes_groups",
        n_genes=N_TOP_MARKER_GENES,
    )

    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    rows = []
    for grp in groups:
        for rank in range(N_TOP_MARKER_GENES):
            rows.append({
                "Cluster": grp,
                "Rank": rank + 1,
                "Gene": result["names"][grp][rank],
                "Score": round(float(result["scores"][grp][rank]), 3),
                "logFC": round(float(result["logfoldchanges"][grp][rank]), 3),
                "padj": float(result["pvals_adj"][grp][rank]),
            })

    marker_df = pd.DataFrame(rows)
    marker_df.to_csv(Path(out_dir) / "cluster_markers.csv", index=False)
    return adata


# ==============================================================================
#  STEP 5 — CONSENSUS TABLE & ANNOTATION
# ==============================================================================

def print_annotation_table(adata):
    """Print a per-cluster summary showing both automated labels."""
    clusters = sorted(
        adata.obs[CLUSTER_KEY].unique().tolist(),
        key=lambda x: int(x) if x.isdigit() else x,
    )

    print(f"\n{'Cluster':<10} {'Cells':>7}  {'CellTypist':<30}  {'Marker Score':<30}")
    print(f"{'-'*10} {'-'*7}  {'-'*30}  {'-'*30}")

    for c in clusters:
        mask = adata.obs[CLUSTER_KEY] == c
        ct = adata.obs.loc[mask, "celltypist_label"].mode()[0] \
            if "celltypist_label" in adata.obs else "n/a"
        ms = adata.obs.loc[mask, "marker_score_label"].mode()[0] \
            if "marker_score_label" in adata.obs else "n/a"
        n = int(mask.sum())
        override = f"  <- MANUAL: {MANUAL_ANNOTATION[str(c)]}" if str(c) in MANUAL_ANNOTATION else ""
        print(f"{str(c):<10} {n:>7,}  {ct:<30}  {ms:<30}{override}")

    print("\nTo override a label, add it to MANUAL_ANNOTATION and re-run.")


def apply_annotation(adata):
    clusters = adata.obs[CLUSTER_KEY].unique().tolist()
    mapping = {}
    for c in clusters:
        if str(c) in MANUAL_ANNOTATION:
            mapping[str(c)] = MANUAL_ANNOTATION[str(c)]
        elif "celltypist_label" in adata.obs:
            mask = adata.obs[CLUSTER_KEY] == c
            mapping[str(c)] = adata.obs.loc[mask, "celltypist_label"].mode()[0]
        else:
            mapping[str(c)] = f"Cluster_{c}"

    adata.obs["cell_type"] = (
        adata.obs[CLUSTER_KEY]
        .astype(str)
        .map(mapping)
        .fillna("Unknown")
        .astype("category")
    )
    return adata


def save_umap_plot(adata, plot_dir):
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    if "X_umap" not in adata.obsm:
        return

    umap = adata.obsm["X_umap"]
    obs = adata.obs.copy()
    obs["UMAP1"] = umap[:, 0]
    obs["UMAP2"] = umap[:, 1]

    fig = px.scatter(
        obs, x="UMAP1", y="UMAP2", color="cell_type",
        title="UMAP — Final Cell Type Annotation",
        template="plotly_white", opacity=0.75,
    )
    fig.update_traces(marker_size=3)
    fig.update_layout(height=540)
    fig.write_html(plot_dir / "01_umap_cell_types.html")


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    out_dir = Path(OUTPUT_PATH).parent
    plot_dir = out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    adata = load_data()
    adata = run_celltypist(adata)
    adata = run_marker_scoring(adata)
    adata = find_cluster_markers(adata, out_dir)
    print_annotation_table(adata)
    adata = apply_annotation(adata)
    save_umap_plot(adata, plot_dir)

    save_adata(adata, OUTPUT_PATH)
    print(f"Saved → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
