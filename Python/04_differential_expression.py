#!/usr/bin/env python3
"""
04_differential_expression.py
==============================
Pseudobulk differential expression and pathway enrichment.

WHAT THIS SCRIPT DOES:
  1. Cell proportion analysis — do cell type abundances shift between conditions?
  2. Pseudobulk aggregation — sum raw counts per (sample × cell type)
  3. PyDESeq2 — negative binomial DE test per cell type
  4. Volcano plot saved as HTML
  5. Pathway enrichment via GSEApy / Enrichr (requires internet)

IMPORTANT — WHY PSEUDOBULK?
  Direct cell-level tests (Wilcoxon) treat each cell as an independent
  observation, inflating significance due to pseudoreplication. Pseudobulk
  aggregates cells by biological replicate (sample) first, which is
  statistically correct. You need >= 3 samples per condition.

HOW TO USE:
  1. Edit CONFIGURATION below.
  2. Run:  python 04_differential_expression.py

OUTPUT:
  results/04_de/de_results.csv
  results/04_de/pathway_enrichment.csv
  results/04_de/plots/
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

# Column names in adata.obs
CELL_TYPE_COL = "cell_type"
SAMPLE_COL = "sample"
CONDITION_COL = "condition"

# --- Comparison ---------------------------------------------------------------
REFERENCE_CONDITION = "Control"    # denominator (e.g., control/untreated)
TEST_CONDITION = "Treatment"       # numerator (e.g., treated)

# Leave empty ("") to run DE for all cell types.
TARGET_CELL_TYPE = ""

# --- DE thresholds ------------------------------------------------------------
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 0.5
MIN_PSEUDOBULK_CELLS = 10

# --- Pathway enrichment -------------------------------------------------------
RUN_ENRICHR = True
ENRICHR_GENE_SETS = ["MSigDB_Hallmark_2020", "KEGG_2021_Human"]

# --- Output ------------------------------------------------------------------
OUTPUT_DIR = "results/04_de"


# ==============================================================================
#  STEP 1 — LOAD & VALIDATE
# ==============================================================================

def load_data():
    p = Path(INPUT_PATH)
    if not p.exists():
        print(f"File not found: {p}")
        sys.exit(1)
    adata = sc.read_h5ad(p)

    missing = [c for c in [CELL_TYPE_COL, SAMPLE_COL, CONDITION_COL]
               if c not in adata.obs.columns]
    if missing:
        print(f"Missing columns: {', '.join(missing)}")
        print(f"Available: {', '.join(adata.obs.columns.tolist())}")
        sys.exit(1)

    return adata


# ==============================================================================
#  STEP 2 — CELL PROPORTION ANALYSIS
# ==============================================================================

def cell_proportion_analysis(adata, plot_dir):
    obs = adata.obs[[CELL_TYPE_COL, SAMPLE_COL, CONDITION_COL]].copy()
    counts = (
        obs.groupby([SAMPLE_COL, CELL_TYPE_COL])
        .size()
        .reset_index(name="n_cells")
    )
    totals = counts.groupby(SAMPLE_COL)["n_cells"].transform("sum")
    counts["proportion"] = counts["n_cells"] / totals
    sample_cond = obs[[SAMPLE_COL, CONDITION_COL]].drop_duplicates().set_index(SAMPLE_COL)
    counts[CONDITION_COL] = counts[SAMPLE_COL].map(sample_cond[CONDITION_COL])

    fig = px.bar(
        counts,
        x=SAMPLE_COL, y="proportion", color=CELL_TYPE_COL,
        barmode="stack", facet_col=CONDITION_COL,
        title="Cell Type Proportions per Sample",
        template="plotly_white",
        labels={"proportion": "Fraction of cells"},
    )
    fig.update_layout(height=480)
    fig.write_html(Path(plot_dir) / "01_cell_proportions.html")


# ==============================================================================
#  STEP 3 — PSEUDOBULK DE (PyDESeq2)
# ==============================================================================

def run_pseudobulk_de(adata):
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("pydeseq2 not installed (pip install pydeseq2).")
        sys.exit(1)

    all_types = adata.obs[CELL_TYPE_COL].unique().tolist()
    types_to_run = [TARGET_CELL_TYPE] if TARGET_CELL_TYPE else all_types

    all_results = []

    for ct in types_to_run:
        mask = adata.obs[CELL_TYPE_COL] == ct
        sub = adata[mask]

        pb_counts = {}
        pb_meta = {}
        for s in sub.obs[SAMPLE_COL].unique():
            s_mask = sub.obs[SAMPLE_COL] == s
            if s_mask.sum() < MIN_PSEUDOBULK_CELLS:
                continue
            X = sub.X[s_mask]
            if hasattr(X, "toarray"):
                X = X.toarray()
            pb_counts[s] = X.sum(axis=0).astype(int)
            cond = sub.obs.loc[s_mask, CONDITION_COL].iloc[0]
            pb_meta[s] = {CONDITION_COL: cond}

        if len(pb_counts) < 2:
            continue

        counts_df = pd.DataFrame(pb_counts, index=adata.var_names).T
        meta_df = pd.DataFrame(pb_meta).T[[CONDITION_COL]]

        keep = meta_df[CONDITION_COL].isin([REFERENCE_CONDITION, TEST_CONDITION])
        counts_df = counts_df[keep]
        meta_df = meta_df[keep]

        if meta_df[CONDITION_COL].nunique() < 2:
            continue

        counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]

        try:
            dds = DeseqDataSet(
                counts=counts_df, metadata=meta_df,
                design_factors=CONDITION_COL, refit_cooks=True, quiet=True,
            )
            dds.deseq2()

            stat = DeseqStats(
                dds,
                contrast=[CONDITION_COL, TEST_CONDITION, REFERENCE_CONDITION],
                quiet=True,
            )
            stat.summary()
            res = stat.results_df.copy()
            res["cell_type"] = ct
            all_results.append(res)
        except Exception as e:
            print(f"{ct}: DE failed — {e}")

    if not all_results:
        print("No DE results. Check REFERENCE_CONDITION/TEST_CONDITION match adata.obs values.")
        sys.exit(1)

    de_results = pd.concat(all_results).reset_index().rename(columns={"index": "gene"})
    de_results["significant"] = (
        (de_results["padj"] <= PADJ_CUTOFF)
        & (de_results["log2FoldChange"].abs() >= LFC_CUTOFF)
    )
    return de_results


# ==============================================================================
#  STEP 4 — PLOTS
# ==============================================================================

def save_volcano_plot(de_results, plot_dir):
    df = de_results.copy()
    df["-log10padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["colour"] = "NS"
    df.loc[(df["log2FoldChange"] > LFC_CUTOFF) & (df["padj"] <= PADJ_CUTOFF), "colour"] = "Up"
    df.loc[(df["log2FoldChange"] < -LFC_CUTOFF) & (df["padj"] <= PADJ_CUTOFF), "colour"] = "Down"

    fig = px.scatter(
        df, x="log2FoldChange", y="-log10padj",
        color="colour",
        color_discrete_map={"Up": "#C44E52", "Down": "#4C72B0", "NS": "#AAAAAA"},
        hover_data=["gene", "cell_type"],
        facet_col="cell_type" if df["cell_type"].nunique() > 1 else None,
        title=f"Volcano Plot: {TEST_CONDITION} vs {REFERENCE_CONDITION}",
        template="plotly_white", opacity=0.6,
        labels={"log2FoldChange": "log₂FC", "-log10padj": "-log₁₀(padj)"},
    )
    fig.update_traces(marker_size=4)
    fig.update_layout(height=480)
    fig.write_html(Path(plot_dir) / "02_volcano.html")


# ==============================================================================
#  STEP 5 — PATHWAY ENRICHMENT
# ==============================================================================

def run_pathway_enrichment(de_results):
    if not RUN_ENRICHR:
        return None

    try:
        import gseapy as gp
    except ImportError:
        print("gseapy not installed (pip install gseapy). Skipping enrichment.")
        return None

    up_genes = de_results.loc[
        (de_results["log2FoldChange"] > LFC_CUTOFF) & (de_results["padj"] <= PADJ_CUTOFF),
        "gene",
    ].tolist()
    dn_genes = de_results.loc[
        (de_results["log2FoldChange"] < -LFC_CUTOFF) & (de_results["padj"] <= PADJ_CUTOFF),
        "gene",
    ].tolist()

    enr_results = []
    for direction, genes in [("Up", up_genes), ("Down", dn_genes)]:
        if len(genes) < 5:
            continue
        for lib in ENRICHR_GENE_SETS:
            try:
                enr = gp.enrichr(
                    gene_list=genes, gene_sets=lib,
                    organism="human", outdir=None, verbose=False,
                )
                r = enr.results.head(10).copy()
                r["direction"] = direction
                r["library"] = lib
                enr_results.append(r)
            except Exception as e:
                print(f"{lib} ({direction}) enrichment failed: {e}")

    if not enr_results:
        return None

    return pd.concat(enr_results)


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    out_dir = Path(OUTPUT_DIR)
    plot_dir = out_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    adata = load_data()
    cell_proportion_analysis(adata, plot_dir)
    de_results = run_pseudobulk_de(adata)
    save_volcano_plot(de_results, plot_dir)

    sig = de_results[de_results["significant"]].sort_values("padj").head(20)
    if not sig.empty:
        print(sig[["gene", "cell_type", "log2FoldChange", "padj"]].to_string(index=False))

    enrichment_df = run_pathway_enrichment(de_results)

    de_path = out_dir / "de_results.csv"
    de_results.to_csv(de_path, index=False)

    if enrichment_df is not None:
        enrichment_df.to_csv(out_dir / "pathway_enrichment.csv", index=False)

    print(f"Saved → {out_dir}/")


if __name__ == "__main__":
    main()
