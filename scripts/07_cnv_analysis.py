#!/usr/bin/env python3
"""
07_cnv_analysis.py
==================
Copy number variation (CNV) analysis using inferCNVpy.

WHAT THIS SCRIPT DOES:
  1. Load annotated tumour + normal data
  2. Annotate genes with chromosome positions (via Ensembl or a custom file)
  3. Run inferCNVpy to infer CNV from expression patterns
  4. Classify cells as aneuploid (tumour) vs. diploid (normal)
  5. Save CNV-annotated AnnData

WHEN TO USE:
  Tumour datasets where you need to distinguish malignant cells from normal
  stromal/immune cells, or to identify tumour subclones.

IMPORTANT CAVEATS:
  - inferCNVpy infers CNV from expression — ~80% concordance with DNA WGS
  - Always interpret in context of cell type annotations
  - A matched normal sample dramatically improves specificity

HOW TO USE:
  1. Edit CONFIGURATION below — especially NORMAL_CELL_TYPES and TUMOUR_CELL_TYPES.
  2. Run:  python 07_cnv_analysis.py

OUTPUT:
  results/07_cnv/adata_cnv.h5ad
  results/07_cnv/plots/
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

# --- Reference (normal) and tumour populations --------------------------------
# Normal cells are used as the baseline for CNV inference.
# These should be non-malignant cells (immune, stromal, etc.).
NORMAL_CELL_TYPES = [
    "T cell",
    "B cell",
    "NK cell",
]

# Cell types you suspect are tumour / malignant.
TUMOUR_CELL_TYPES = [
    "Epithelial",
]

# --- Gene position annotation -------------------------------------------------
# Leave empty ("") to auto-download positions from Ensembl via pybiomart.
# Or provide a path to a custom TSV file with columns:
#   gene_name  chromosome  start  end    (no header, tab-separated)
GENE_ORDER_FILE = ""

# Reference genome version for Ensembl lookup.
# Options: "GRCh38", "GRCh37", "GRCm39", "GRCm38"
GENOME_VERSION = "GRCh38"

# --- inferCNVpy parameters ----------------------------------------------------
# Number of genes per sliding window (smoothing bandwidth).
WINDOW_SIZE = 100

# Step size between windows.
STEP_SIZE = 10

# CNV score threshold for classifying a cell as aneuploid.
# Cells with CNV score > (median_normal_score + threshold) are called aneuploid.
CNV_SCORE_THRESHOLD = 0.1

# --- Output ------------------------------------------------------------------
OUTPUT_PATH = "results/07_cnv/adata_cnv.h5ad"


# ==============================================================================
#  STEP 1 — LOAD DATA
# ==============================================================================

def load_data():
    p = Path(INPUT_PATH)
    if not p.exists():
        print(f"File not found: {p}")
        sys.exit(1)
    adata = sc.read_h5ad(p)

    all_types = adata.obs[CELL_TYPE_COL].unique().tolist()
    missing_normal = [t for t in NORMAL_CELL_TYPES if t not in all_types]
    missing_tumour = [t for t in TUMOUR_CELL_TYPES if t not in all_types]

    if missing_normal:
        print(f"Normal types not found: {missing_normal}\nAvailable: {all_types}")
        sys.exit(1)
    if missing_tumour:
        print(f"Tumour types not found: {missing_tumour}\nAvailable: {all_types}")
        sys.exit(1)

    return adata


# ==============================================================================
#  STEP 2 — GENE CHROMOSOME POSITIONS
# ==============================================================================

def annotate_gene_positions(adata):
    if GENE_ORDER_FILE:
        gf = Path(GENE_ORDER_FILE)
        if not gf.exists():
            print(f"Gene order file not found: {gf}")
            sys.exit(1)
        gene_pos = pd.read_csv(gf, sep="\t", header=None,
                               names=["gene_name", "chromosome", "start", "end"],
                               index_col=0)
    else:
        try:
            import pybiomart as pbm
            dataset_name = {
                "GRCh38": "hsapiens_gene_ensembl",
                "GRCh37": "hsapiens_gene_ensembl",
                "GRCm39": "mmusculus_gene_ensembl",
                "GRCm38": "mmusculus_gene_ensembl",
            }[GENOME_VERSION]
            server = pbm.Server(host="http://www.ensembl.org")
            mart = server["ENSEMBL_MART_ENSEMBL"][dataset_name]
            gene_pos = mart.query(
                attributes=["hgnc_symbol", "chromosome_name", "start_position", "end_position"]
            ).dropna()
            gene_pos.columns = ["gene_name", "chromosome", "start", "end"]
            gene_pos = (
                gene_pos[gene_pos["gene_name"].isin(adata.var_names)]
                .drop_duplicates("gene_name")
                .set_index("gene_name")
            )
        except ImportError:
            print("pybiomart not installed (pip install pybiomart). "
                  "Provide a custom GENE_ORDER_FILE instead.")
            sys.exit(1)
        except Exception as e:
            print(f"Ensembl download failed: {e}\n"
                  "Check internet connection, or provide a GENE_ORDER_FILE.")
            sys.exit(1)

    for col in ["chromosome", "start", "end"]:
        if col in gene_pos.columns:
            adata.var[col] = gene_pos.reindex(adata.var_names)[col]

    return adata


# ==============================================================================
#  STEP 3 — RUN inferCNVpy
# ==============================================================================

def run_infercnv(adata):
    try:
        import infercnvpy as cnv
    except ImportError:
        print("infercnvpy not installed (pip install infercnvpy).")
        sys.exit(1)

    all_interest = NORMAL_CELL_TYPES + TUMOUR_CELL_TYPES
    mask = adata.obs[CELL_TYPE_COL].isin(all_interest)
    adata_cnv = adata[mask].copy()

    REF_KEY = "is_reference"
    adata_cnv.obs[REF_KEY] = adata_cnv.obs[CELL_TYPE_COL].isin(NORMAL_CELL_TYPES)

    if "counts" in adata_cnv.layers:
        import scipy.sparse as sp
        adata_cnv.X = adata_cnv.layers["counts"]

    cnv.tl.infercnv(
        adata_cnv,
        reference_key=REF_KEY,
        reference_cat=True,
        window_size=WINDOW_SIZE,
        step=STEP_SIZE,
        key_added="cnv",
    )
    cnv.tl.pca(adata_cnv, use_rep="X_cnv")
    cnv.tl.umap(adata_cnv)
    cnv.tl.leiden(adata_cnv)

    cnv_mat = adata_cnv.obsm["X_cnv"]
    adata_cnv.obs["cnv_score"] = np.var(
        cnv_mat.toarray() if hasattr(cnv_mat, "toarray") else cnv_mat,
        axis=1,
    )

    med_ref = adata_cnv.obs.loc[adata_cnv.obs[REF_KEY], "cnv_score"].median()
    adata_cnv.obs["cnv_status"] = np.where(
        adata_cnv.obs["cnv_score"] > med_ref + CNV_SCORE_THRESHOLD,
        "Aneuploid",
        "Diploid",
    )

    adata.obs.loc[adata_cnv.obs.index, "cnv_score"] = adata_cnv.obs["cnv_score"]
    adata.obs.loc[adata_cnv.obs.index, "cnv_status"] = adata_cnv.obs["cnv_status"]

    n_aneupl = (adata_cnv.obs["cnv_status"] == "Aneuploid").sum()
    pct = 100 * n_aneupl / len(adata_cnv.obs)
    print(f"Aneuploid cells: {n_aneupl:,} ({pct:.1f}%)")
    return adata, adata_cnv


# ==============================================================================
#  STEP 4 — VISUALISATION
# ==============================================================================

def save_cnv_plots(adata, adata_cnv, plot_dir):
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    if "X_umap" in adata.obsm and "cnv_status" in adata.obs.columns:
        umap = adata.obsm["X_umap"]
        obs = adata.obs.copy()
        obs["UMAP1"] = umap[:, 0]
        obs["UMAP2"] = umap[:, 1]
        fig = px.scatter(
            obs, x="UMAP1", y="UMAP2", color="cnv_status",
            color_discrete_map={"Aneuploid": "#C44E52", "Diploid": "#4C72B0"},
            title="CNV Status on Original UMAP",
            template="plotly_white", opacity=0.7,
        )
        fig.update_traces(marker_size=3)
        fig.update_layout(height=480)
        fig.write_html(plot_dir / "01_umap_cnv_status.html")

    try:
        import infercnvpy as cnv
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(1, 2, figsize=(14, 5))
        cnv.pl.umap(adata_cnv, color="cnv_status", ax=axs[0], show=False,
                    title="CNV Status")
        cnv.pl.umap(adata_cnv, color="cnv_score", ax=axs[1], show=False,
                    title="CNV Score", color_map="Reds")
        plt.tight_layout()
        fig.savefig(plot_dir / "02_cnv_umap.png", bbox_inches="tight", dpi=110)
        plt.close(fig)

        fig2 = cnv.pl.chromosome_heatmap(adata_cnv, groupby="cnv_status",
                                          return_fig=True)
        fig2.savefig(plot_dir / "03_chromosome_heatmap.png", bbox_inches="tight", dpi=100)
        plt.close(fig2)
    except Exception as e:
        print(f"CNV plot generation failed: {e}")


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    plot_dir = Path(OUTPUT_PATH).parent / "plots"

    adata = load_data()
    adata = annotate_gene_positions(adata)
    adata, adata_cnv = run_infercnv(adata)
    save_cnv_plots(adata, adata_cnv, plot_dir)

    save_adata(adata_cnv, OUTPUT_PATH)
    print(f"Saved → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
