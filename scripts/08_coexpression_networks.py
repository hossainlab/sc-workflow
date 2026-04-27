#!/usr/bin/env python3
"""
08_coexpression_networks.py
============================
Gene co-expression network analysis using PyWGCNA.

WHAT THIS SCRIPT DOES:
  1. Subset to target cell type and build metacells (k-NN aggregation)
     — metacells average k nearest neighbours to fill technical zeros
  2. Soft-thresholding power selection (scale-free topology test)
  3. WGCNA network construction and module detection
  4. Hub gene identification per module (by intramodular connectivity)
  5. Pathway enrichment per module via Enrichr
  6. Save all results

WHY METACELLS?
  Single-cell data is 80-90% sparse (zeros). Averaging k nearest neighbours
  in PCA space fills in technical zeros without losing cell-type specificity.

HOW TO USE:
  1. Edit CONFIGURATION below.
  2. Run:  python 08_coexpression_networks.py

OUTPUT:
  results/08_coexpression/module_assignments.csv
  results/08_coexpression/module_eigengenes.csv
  results/08_coexpression/hub_genes.csv
  results/08_coexpression/module_enrichment.csv
  results/08_coexpression/plots/
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
SAMPLE_COL = "sample"
CONDITION_COL = "condition"   # used for module-trait correlation

# --- Target cell type ---------------------------------------------------------
# Which cell type to build the network for.
# Leave empty ("") to use all cells (not recommended — very slow and noisy).
TARGET_CELL_TYPE = "CD4+ T cell"

# --- Metacell parameters ------------------------------------------------------
# k = number of nearest neighbours to average per metacell.
# Rule of thumb: k = 20-30 for most datasets.
METACELL_K = 20

# Two metacells sharing more than this many cells are not both kept.
MAX_SHARED = 10

# Only keep genes expressed (count > 0) in at least this fraction of metacells.
MIN_EXPR_FRACTION = 0.05

# --- WGCNA parameters --------------------------------------------------------
# Soft-thresholding power β.
# Set to 0 to auto-detect (chooses lowest β where scale-free R² ≥ 0.85).
SOFT_POWER = 0

# Minimum number of genes per module.
MIN_MODULE_SIZE = 30

# Modules with eigengene correlation > (1 - merge_cut_height) are merged.
MERGE_CUT_HEIGHT = 0.25

# --- Pathway enrichment per module --------------------------------------------
RUN_MODULE_ENRICHMENT = True
MODULE_TOP_N_GENES = 50   # top genes per module to submit to Enrichr

# --- Output ------------------------------------------------------------------
OUTPUT_DIR = "results/08_coexpression"


# ==============================================================================
#  HELPERS
# ==============================================================================

def print_header(title):
    print()
    print("=" * 60)
    print(f"  {title}")
    print("=" * 60)


def print_ok(msg):
    print(f"  [OK]  {msg}")


def print_info(msg):
    print(f"  [--]  {msg}")


# ==============================================================================
#  STEP 1 — LOAD & SUBSET
# ==============================================================================

def load_and_subset():
    print_header("Step 1 — Load & Subset Data")
    p = Path(INPUT_PATH)
    if not p.exists():
        print(f"  [!!]  File not found: {p}")
        sys.exit(1)
    adata = sc.read_h5ad(p)
    print_ok(f"Loaded {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if TARGET_CELL_TYPE:
        mask = adata.obs[CELL_TYPE_COL] == TARGET_CELL_TYPE
        adata_sub = adata[mask].copy()
        if adata_sub.n_obs < 100:
            print(f"  [!!]  Only {adata_sub.n_obs} cells for '{TARGET_CELL_TYPE}'.")
            print("        Need >= 100 cells (>=200 recommended).")
            sys.exit(1)
        print_ok(f"Subset to '{TARGET_CELL_TYPE}': {adata_sub.n_obs:,} cells")
    else:
        adata_sub = adata.copy()
        print_info("Using all cells (TARGET_CELL_TYPE is empty)")

    return adata_sub


# ==============================================================================
#  STEP 2 — METACELL CONSTRUCTION
# ==============================================================================

def build_metacells(adata_sub):
    print_header("Step 2 — Metacell Construction")
    print_info(f"k={METACELL_K} nearest neighbours per metacell")

    from sklearn.neighbors import NearestNeighbors
    import scipy.sparse as sp

    # Use raw counts if available
    X = adata_sub.layers["counts"] if "counts" in adata_sub.layers else adata_sub.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(float)

    samples = (adata_sub.obs[SAMPLE_COL].unique().tolist()
               if SAMPLE_COL in adata_sub.obs.columns else ["all"])

    metacell_mats = []
    metacell_meta = []

    for s in samples:
        if s == "all":
            s_mask = np.ones(adata_sub.n_obs, dtype=bool)
        else:
            s_mask = (adata_sub.obs[SAMPLE_COL] == s).values

        X_s = X[s_mask]
        pca_s = (adata_sub.obsm["X_pca"][s_mask]
                 if "X_pca" in adata_sub.obsm else X_s[:, :50])

        n_mc_target = max(1, X_s.shape[0] // METACELL_K)
        nbrs = NearestNeighbors(n_neighbors=METACELL_K, algorithm="ball_tree").fit(pca_s)
        _, indices = nbrs.kneighbors(pca_s)

        used = set()
        for i in range(X_s.shape[0]):
            knn = indices[i]
            overlap = sum(j in used for j in knn)
            if overlap > MAX_SHARED:
                continue
            metacell_mats.append(X_s[knn].mean(axis=0))
            metacell_meta.append({"sample": s})
            used.update(knn.tolist())
            if len(metacell_mats) >= n_mc_target:
                break

    if not metacell_mats:
        print("  [!!]  Metacell construction failed. Try reducing METACELL_K.")
        sys.exit(1)

    mc_matrix = np.vstack(metacell_mats)   # metacells × genes
    mc_meta = pd.DataFrame(metacell_meta)

    # Gene filter
    expr_frac = (mc_matrix > 0).mean(axis=0)
    gene_mask = expr_frac >= MIN_EXPR_FRACTION
    mc_matrix = mc_matrix[:, gene_mask]
    genes_kept = adata_sub.var_names[gene_mask].tolist()

    metacell_df = pd.DataFrame(mc_matrix, columns=genes_kept)
    print_ok(f"Built {metacell_df.shape[0]:,} metacells × {metacell_df.shape[1]:,} genes")
    print_ok(f"(Genes expressed in >= {MIN_EXPR_FRACTION*100:.0f}% of metacells retained)")
    return metacell_df, mc_meta, genes_kept


# ==============================================================================
#  STEP 3 — SOFT-THRESHOLDING POWER SELECTION
# ==============================================================================

def select_soft_power(metacell_df):
    print_header("Step 3 — Soft-Thresholding Power Selection")

    if SOFT_POWER > 0:
        print_ok(f"Using manually set SOFT_POWER = {SOFT_POWER}")
        return SOFT_POWER

    powers = list(range(1, 31))
    r2_scores = []

    corr = np.corrcoef(metacell_df.values.T)   # genes × genes
    np.fill_diagonal(corr, 0)
    corr_abs = np.abs(corr)

    for beta in powers:
        adj = corr_abs ** beta
        k = adj.sum(axis=1) - 1   # connectivity degree
        log_k = np.log10(k + 1e-10)
        hist_counts, bin_edges = np.histogram(k, bins=20)
        log_pk = np.log10(hist_counts / (len(k) + 1e-10) + 1e-10)
        x = np.linspace(log_k.min(), log_k.max(), 20)
        valid = np.isfinite(log_pk) & (log_pk > -np.inf)
        if valid.sum() > 3:
            fit = np.polyfit(x[valid], log_pk[valid], 1)
            yhat = np.polyval(fit, x[valid])
            ss_res = ((log_pk[valid] - yhat) ** 2).sum()
            ss_tot = ((log_pk[valid] - log_pk[valid].mean()) ** 2).sum()
            r2 = max(0.0, 1 - ss_res / (ss_tot + 1e-10))
        else:
            r2 = 0.0
        r2_scores.append(r2)

    sft_df = pd.DataFrame({"Power": powers, "R2": r2_scores})
    auto_candidates = sft_df[sft_df["R2"] >= 0.85]
    auto_power = int(auto_candidates["Power"].iloc[0]) if len(auto_candidates) > 0 else 6

    print()
    print("  Power  R²")
    print("  ─────  ────")
    for _, row in sft_df.iterrows():
        marker = " <-- selected" if row["Power"] == auto_power else ""
        print(f"  {int(row['Power']):>5}  {row['R2']:.3f}{marker}")

    print_ok(f"Auto-selected soft power β = {auto_power} "
             f"(R² = {sft_df.loc[sft_df['Power']==auto_power, 'R2'].values[0]:.3f})")
    return auto_power


# ==============================================================================
#  STEP 4 — WGCNA MODULE DETECTION
# ==============================================================================

def run_wgcna(metacell_df, mc_meta, selected_power):
    print_header("Step 4 — WGCNA Network Construction & Module Detection")

    try:
        import PyWGCNA as pw
    except ImportError:
        print("  [!!]  PyWGCNA not installed. Run: pip install PyWGCNA")
        sys.exit(1)

    name = (TARGET_CELL_TYPE.replace(" ", "_").replace("+", "pos")
            if TARGET_CELL_TYPE else "all_cells")

    wgcna = pw.WGCNA(
        name=name,
        species="human",
        geneExp=metacell_df.T,   # PyWGCNA expects genes × samples
        TPMcutoff=0,
        save=False,
    )

    wgcna.preprocess()
    wgcna.findModules(
        power=selected_power,
        minModuleSize=MIN_MODULE_SIZE,
        deepSplit=4,
        mergeCutHeight=MERGE_CUT_HEIGHT,
    )
    wgcna.sampleInfo = mc_meta.copy()

    n_modules = len(wgcna.moduleColors.unique()) - 1  # exclude grey (unassigned)
    print_ok(f"WGCNA complete: {n_modules} co-expression modules detected "
             f"across {metacell_df.shape[1]:,} genes")
    return wgcna


# ==============================================================================
#  STEP 5 — HUB GENES & VISUALISATION
# ==============================================================================

def get_hub_genes(wgcna):
    print_header("Step 5 — Hub Genes (Intramodular Connectivity)")

    kme = wgcna.kMEtable
    hub_genes = (
        kme.reset_index()
        .melt(id_vars="index", var_name="module", value_name="kME")
        .rename(columns={"index": "gene"})
        .sort_values("kME", ascending=False)
        .groupby("module")
        .head(10)
        .sort_values(["module", "kME"], ascending=[True, False])
    )

    print()
    print("  Top 3 hub genes per module:")
    for module, group in hub_genes.groupby("module"):
        top3 = group.head(3)["gene"].tolist()
        print(f"    {module}: {', '.join(top3)}")

    return hub_genes


def save_eigengene_plot(wgcna, plot_dir):
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    me = wgcna.MEs
    fig = px.imshow(
        me.corr(),
        color_continuous_scale="RdBu_r",
        zmin=-1, zmax=1,
        title="Module Eigengene Correlation Matrix",
        template="plotly_white",
    )
    fig.update_layout(height=500)
    fig.write_html(plot_dir / "01_eigengene_correlation.html")
    print_ok(f"Eigengene plot saved → {plot_dir}/01_eigengene_correlation.html")


# ==============================================================================
#  STEP 6 — PATHWAY ENRICHMENT PER MODULE
# ==============================================================================

def run_module_enrichment(wgcna):
    print_header("Step 6 — Pathway Enrichment per Module (Enrichr)")

    if not RUN_MODULE_ENRICHMENT:
        print_info("Module enrichment skipped (RUN_MODULE_ENRICHMENT = False)")
        return None

    try:
        import gseapy as gp
    except ImportError:
        print("  [!!]  gseapy not installed. Run: pip install gseapy")
        return None

    results = []
    modules = [m for m in wgcna.moduleColors.unique() if m != "grey"]

    for mod in modules:
        genes = wgcna.datExpr.columns[wgcna.moduleColors == mod].tolist()
        genes = genes[:MODULE_TOP_N_GENES]
        if len(genes) < 5:
            continue
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=["MSigDB_Hallmark_2020", "KEGG_2021_Human"],
                organism="human",
                outdir=None,
                verbose=False,
            )
            r = enr.results.head(5).copy()
            r["module"] = mod
            r["n_genes"] = len(genes)
            results.append(r)
            print_ok(f"Module {mod}: {len(r)} enrichment terms")
        except Exception as e:
            print_info(f"Module {mod}: enrichment failed — {e}")

    if not results:
        print_info("No enrichment results returned")
        return None

    return pd.concat(results)


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    print()
    print("=" * 60)
    print("  scRNA-seq Co-expression Networks (PyWGCNA)")
    print("=" * 60)
    print(f"  Input       : {INPUT_PATH}")
    print(f"  Cell type   : {TARGET_CELL_TYPE or 'all cells'}")
    print(f"  Metacell k  : {METACELL_K}")
    print(f"  Output      : {OUTPUT_DIR}")

    out_dir = Path(OUTPUT_DIR)
    plot_dir = out_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    adata_sub = load_and_subset()
    metacell_df, mc_meta, genes_kept = build_metacells(adata_sub)
    selected_power = select_soft_power(metacell_df)
    wgcna = run_wgcna(metacell_df, mc_meta, selected_power)
    hub_genes = get_hub_genes(wgcna)
    save_eigengene_plot(wgcna, plot_dir)
    enrichment_df = run_module_enrichment(wgcna)

    print_header("Save")

    # Module assignments
    module_path = out_dir / "module_assignments.csv"
    pd.DataFrame({
        "gene": wgcna.datExpr.columns,
        "module": wgcna.moduleColors,
    }).to_csv(module_path, index=False)
    print_ok(f"Module assignments → {module_path}")

    # Module eigengenes
    me_path = out_dir / "module_eigengenes.csv"
    wgcna.MEs.to_csv(me_path)
    print_ok(f"Module eigengenes → {me_path}")

    # Hub genes
    hub_path = out_dir / "hub_genes.csv"
    hub_genes.to_csv(hub_path, index=False)
    print_ok(f"Hub genes → {hub_path}")

    # Pathway enrichment
    if enrichment_df is not None:
        enr_path = out_dir / "module_enrichment.csv"
        enrichment_df.to_csv(enr_path, index=False)
        print_ok(f"Module enrichment → {enr_path}")

    print()
    print("  Pipeline complete. All results saved to:", out_dir)
    print()


if __name__ == "__main__":
    main()
