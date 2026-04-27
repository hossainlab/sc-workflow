#!/usr/bin/env python3
"""
06_cell_communication.py
=========================
Infer intercellular signalling using LIANA (aggregates multiple methods).

WHAT THIS SCRIPT DOES:
  1. Run LIANA aggregate score (combines CellChat, NATMI, logFC, and others)
  2. Save top interaction results and heatmap
  3. Optional: run LIANA per condition for comparative analysis
  4. Save results as CSV

WHY LIANA?
  It runs multiple ligand-receptor scoring methods simultaneously and
  aggregates their rankings, giving more robust results than any single method.

IMPORTANT CAVEAT:
  LIANA predicts *potential* communication from co-expression of ligands and
  receptors. Results are hypotheses that require wet-lab validation.

HOW TO USE:
  1. Edit CONFIGURATION below.
  2. Run:  python 06_cell_communication.py

OUTPUT:
  results/06_communication/liana_interactions.csv
  results/06_communication/liana_comparative.csv  (if RUN_COMPARATIVE = True)
  results/06_communication/plots/
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px

sys.path.insert(0, str(Path("..") / "src"))
from sc_workflow.utils import setup_scanpy, save_adata

warnings.filterwarnings("ignore")
setup_scanpy()


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH = "results/03_annotation/adata_annotated.h5ad"
CELL_TYPE_COL = "cell_type"

# Column for condition comparison — leave empty ("") for single-condition analysis.
CONDITION_COL = "condition"

# --- LIANA parameters --------------------------------------------------------
# Methods to aggregate. Options: "cellchat", "connectome", "logfc", "natmi", "singlecellsignalr"
LIANA_METHODS = ["cellchat", "connectome", "logfc", "natmi"]

# Ligand-receptor resource. Options: "consensus", "CellChatDB", "CellPhoneDBv5", "NATMI_LR"
LIANA_RESOURCE = "consensus"

# Minimum fraction of cells in a group that must express the ligand/receptor.
MIN_EXPR_PROP = 0.1

# Number of top interactions to display in the heatmap.
N_TOP_INTERACTIONS = 30

# --- Comparative analysis ----------------------------------------------------
# Set True to run LIANA separately per condition (requires CONDITION_COL to be set).
RUN_COMPARATIVE = False

# --- Output ------------------------------------------------------------------
OUTPUT_DIR = "results/06_communication"


# ==============================================================================
#  STEP 1 — LOAD DATA
# ==============================================================================

def load_data():
    p = Path(INPUT_PATH)
    if not p.exists():
        print(f"File not found: {p}")
        sys.exit(1)
    adata = sc.read_h5ad(p)

    if CELL_TYPE_COL not in adata.obs.columns:
        print(f"Column '{CELL_TYPE_COL}' not found in adata.obs.")
        sys.exit(1)

    return adata


# ==============================================================================
#  STEP 2 — LIANA AGGREGATE SCORE
# ==============================================================================

def run_liana(adata):
    try:
        import liana as li
    except ImportError:
        print("liana-py not installed (pip install liana-py).")
        sys.exit(1)

    adata_li = adata.copy()
    sc.pp.normalize_total(adata_li, target_sum=1e4)
    sc.pp.log1p(adata_li)

    li.mt.rank_aggregate(
        adata_li,
        groupby=CELL_TYPE_COL,
        resource_name=LIANA_RESOURCE,
        methods=LIANA_METHODS,
        expr_prop=MIN_EXPR_PROP,
        verbose=True,
        use_raw=False,
        inplace=True,
    )

    adata.uns["liana_res"] = adata_li.uns["liana_res"]
    return adata


# ==============================================================================
#  STEP 3 — VISUALISATION
# ==============================================================================

def save_liana_plots(adata, plot_dir):
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    if "liana_res" not in adata.uns:
        return

    res = (
        adata.uns["liana_res"]
        .sort_values("aggregate_rank")
        .head(N_TOP_INTERACTIONS)
    )

    heatmap = (
        res.groupby(["source", "target"])["aggregate_rank"]
        .min()
        .reset_index()
        .pivot(index="source", columns="target", values="aggregate_rank")
        .fillna(1.0)
    )
    fig = px.imshow(
        heatmap,
        color_continuous_scale="RdYlBu_r",
        title=f"Top {N_TOP_INTERACTIONS} Interactions: Communication Strength\n"
              "(lower aggregate rank = stronger signal)",
        template="plotly_white",
        text_auto=".2f",
    )
    fig.update_layout(height=500)
    fig.write_html(plot_dir / "01_communication_heatmap.html")

    print("\nTop 20 interactions:")
    print(
        res[["source", "target", "ligand_complex", "receptor_complex", "aggregate_rank"]]
        .head(20)
        .to_string(index=False)
    )


# ==============================================================================
#  STEP 4 — COMPARATIVE ANALYSIS (OPTIONAL)
# ==============================================================================

def run_comparative_liana(adata):
    if not RUN_COMPARATIVE:
        return None

    if not CONDITION_COL or CONDITION_COL not in adata.obs.columns:
        print(f"Condition column '{CONDITION_COL}' not found in adata.obs.")
        return None

    try:
        import liana as li
    except ImportError:
        print("liana-py not installed.")
        return None

    conditions = adata.obs[CONDITION_COL].unique().tolist()
    results = {}

    for cond in conditions:
        mask = adata.obs[CONDITION_COL] == cond
        sub = adata[mask].copy()
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        li.mt.rank_aggregate(
            sub, groupby=CELL_TYPE_COL, resource_name="consensus",
            expr_prop=MIN_EXPR_PROP, verbose=False, use_raw=False, inplace=True,
        )
        results[cond] = sub.uns["liana_res"].copy()
        results[cond][CONDITION_COL] = cond

    return pd.concat(list(results.values()), ignore_index=True)


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    out_dir = Path(OUTPUT_DIR)
    plot_dir = out_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = load_data()
    adata = run_liana(adata)
    save_liana_plots(adata, plot_dir)
    comparative_df = run_comparative_liana(adata)

    if "liana_res" in adata.uns:
        adata.uns["liana_res"].to_csv(out_dir / "liana_interactions.csv", index=False)

    if comparative_df is not None:
        comparative_df.to_csv(out_dir / "liana_comparative.csv", index=False)

    print(f"Saved → {out_dir}/")


if __name__ == "__main__":
    main()
