#!/usr/bin/env python3
"""
05_trajectory_analysis.py
==========================
Trajectory and pseudotime analysis using CellRank.

WHAT THIS SCRIPT DOES:
  1. Load annotated data and subset to selected cell types / lineage
  2. Compute RNA velocity (scVelo) or diffusion pseudotime as the trajectory driver
  3. CellRank kernel construction
  4. Macrostate identification (initial, terminal, intermediate states)
  5. Fate probability computation for each cell
  6. Save trajectory AnnData

WHEN TO USE:
  Developmental systems (stem cell differentiation), time-course experiments
  (treatment response), or any biology with cell state transitions.

HOW TO USE:
  1. Edit CONFIGURATION below.
  2. Run:  python 05_trajectory_analysis.py

OUTPUT:
  results/05_trajectory/adata_trajectory.h5ad
  results/05_trajectory/plots/
"""

import sys
import random
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
random.seed(42)
np.random.seed(42)


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH = "results/03_annotation/adata_annotated.h5ad"
CELL_TYPE_COL = "cell_type"

# --- Lineage selection --------------------------------------------------------
# Cell types to include in the trajectory.
# Leave empty ([]) to use all cell types.
SUBSET_CELL_TYPES = []
# Example: ["Stem cell", "Progenitor", "Effector T cell"]

# Starting cell type for the trajectory.
# Leave empty ("") for automatic detection.
ROOT_CELL_TYPE = ""

# --- Trajectory method --------------------------------------------------------
# "diffusion_pseudotime" : faster, no RNA velocity data needed (recommended to start)
# "scvelo"               : RNA velocity, requires unspliced/spliced count layers
TRAJECTORY_METHOD = "diffusion_pseudotime"

# --- CellRank parameters ------------------------------------------------------
N_MACROSTATES = 4       # number of end states (macrostates) to identify
N_JOBS = 4              # parallel jobs for CellRank computations
DOWNSAMPLE_FRAC = 0.3   # fraction of cells to use (reduces compute time)
                         # set to 1.0 to use all cells

# --- Output ------------------------------------------------------------------
OUTPUT_PATH = "results/05_trajectory/adata_trajectory.h5ad"


# ==============================================================================
#  STEP 1 — LOAD & SUBSET
# ==============================================================================

def load_and_subset():
    p = Path(INPUT_PATH)
    if not p.exists():
        print(f"File not found: {p}")
        sys.exit(1)
    adata = sc.read_h5ad(p)

    if SUBSET_CELL_TYPES:
        mask = adata.obs[CELL_TYPE_COL].isin(SUBSET_CELL_TYPES)
        adata = adata[mask].copy()

    if DOWNSAMPLE_FRAC < 1.0:
        keep_idx = []
        for ct, grp in adata.obs.groupby(CELL_TYPE_COL):
            n = max(50, int(len(grp) * DOWNSAMPLE_FRAC))
            n = min(n, len(grp))
            keep_idx.extend(np.random.choice(grp.index, size=n, replace=False).tolist())
        adata = adata[adata.obs.index.isin(keep_idx)].copy()

    return adata


# ==============================================================================
#  STEP 2 — TRAJECTORY DRIVER (PSEUDOTIME OR VELOCITY)
# ==============================================================================

def compute_trajectory_driver(adata):
    if TRAJECTORY_METHOD == "diffusion_pseudotime":
        sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_pca", random_state=42)
        sc.tl.diffmap(adata)

        root_idx = adata.obsm["X_diffmap"][:, 1].argmin()

        if ROOT_CELL_TYPE:
            type_mask = adata.obs[CELL_TYPE_COL] == ROOT_CELL_TYPE
            if type_mask.any():
                sub_idx = np.where(type_mask)[0]
                root_idx = sub_idx[
                    adata.obsm["X_diffmap"][sub_idx, 1].argmin()
                ]

        adata.uns["iroot"] = int(root_idx)
        sc.tl.dpt(adata)

    elif TRAJECTORY_METHOD == "scvelo":
        try:
            import scvelo as scv
            scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
            scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
            scv.tl.velocity(adata)
            scv.tl.velocity_graph(adata)
        except ImportError:
            print("scvelo not installed (pip install scvelo). Falling back to diffusion pseudotime.")
            return compute_trajectory_driver_dpt(adata)

    return adata


def compute_trajectory_driver_dpt(adata):
    """Fallback: diffusion pseudotime."""
    sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_pca", random_state=42)
    sc.tl.diffmap(adata)
    root_idx = adata.obsm["X_diffmap"][:, 1].argmin()
    adata.uns["iroot"] = int(root_idx)
    sc.tl.dpt(adata)
    return adata


# ==============================================================================
#  STEP 3 — CELLRANK
# ==============================================================================

def run_cellrank(adata):
    try:
        import cellrank as cr
        from cellrank.kernels import ConnectivityKernel, VelocityKernel, PseudotimeKernel
    except ImportError:
        print("cellrank not installed (pip install 'cellrank[knn]').")
        sys.exit(1)

    if TRAJECTORY_METHOD == "scvelo" and "velocity" in adata.layers:
        vk = VelocityKernel(adata).compute_transition_matrix()
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        kernel = 0.8 * vk + 0.2 * ck
    else:
        ptk = PseudotimeKernel(adata, time_key="dpt_pseudotime")
        ptk.compute_transition_matrix(threshold_scheme="soft", nu=0.5)
        ck = ConnectivityKernel(adata).compute_transition_matrix()
        kernel = 0.8 * ptk + 0.2 * ck

    g = cr.estimators.GPCCA(kernel)
    g.compute_schur(n_components=min(N_MACROSTATES + 2, 20))
    g.compute_macrostates(n_states=N_MACROSTATES, cluster_key=CELL_TYPE_COL)
    g.predict_terminal_states()
    g.compute_fate_probabilities()

    adata = g.adata
    adata.uns["_cellrank_estimator"] = g
    return adata


# ==============================================================================
#  STEP 4 — VISUALISATION
# ==============================================================================

def save_trajectory_plots(adata, plot_dir):
    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    if "X_umap" not in adata.obsm:
        return

    umap = adata.obsm["X_umap"]
    obs = adata.obs.copy()
    obs["UMAP1"] = umap[:, 0]
    obs["UMAP2"] = umap[:, 1]

    if "dpt_pseudotime" in obs.columns:
        fig = px.scatter(
            obs, x="UMAP1", y="UMAP2", color="dpt_pseudotime",
            color_continuous_scale="Viridis",
            title="UMAP — Diffusion Pseudotime",
            template="plotly_white", opacity=0.75,
        )
        fig.update_traces(marker_size=4)
        fig.update_layout(height=460)
        fig.write_html(plot_dir / "01_umap_pseudotime.html")

    fig = px.scatter(
        obs, x="UMAP1", y="UMAP2", color=CELL_TYPE_COL,
        title="UMAP — Cell Types",
        template="plotly_white", opacity=0.75,
    )
    fig.update_traces(marker_size=4)
    fig.update_layout(height=460)
    fig.write_html(plot_dir / "02_umap_celltypes.html")


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    plot_dir = Path(OUTPUT_PATH).parent / "plots"

    adata = load_and_subset()
    adata = compute_trajectory_driver(adata)
    adata = run_cellrank(adata)
    save_trajectory_plots(adata, plot_dir)

    # Remove non-serialisable estimator object before saving
    adata_save = adata.copy()
    adata_save.uns.pop("_cellrank_estimator", None)
    save_adata(adata_save, OUTPUT_PATH)
    print(f"Saved → {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
