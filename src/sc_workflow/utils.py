"""
Shared utilities for the single-cell RNA-seq workflow.

All functions operate on AnnData objects and follow scanpy conventions.
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

logger = logging.getLogger(__name__)


# ── Environment ───────────────────────────────────────────────────────────────

def setup_scanpy(verbosity: int = 2, n_jobs: int = -1, dpi: int = 100) -> None:
    """Configure scanpy global settings."""
    sc.settings.verbosity = verbosity
    sc.settings.n_jobs = n_jobs
    sc.settings.set_figure_params(dpi=dpi, facecolor="white", frameon=False)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")


# ── I/O ───────────────────────────────────────────────────────────────────────

def load_data(
    path: Union[str, Path],
    var_names: str = "gene_symbols",
    sample_name: Optional[str] = None,
) -> sc.AnnData:
    """
    Load single-cell data from multiple formats.

    Supported formats:
    - 10x Genomics Cell Ranger output directory
      (contains filtered_feature_bc_matrix/ or raw_feature_bc_matrix/)
    - 10x MTX directory directly (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
    - 10x HDF5 file (.h5)
    - AnnData HDF5 file (.h5ad)

    Parameters
    ----------
    path : str | Path
        Path to data (file or directory).
    var_names : str
        Use "gene_symbols" or "gene_ids" for 10x data.
    sample_name : str, optional
        If provided, adds 'sample' column to adata.obs.

    Returns
    -------
    sc.AnnData with unique obs/var names.
    """
    path = Path(path)

    if path.suffix == ".h5ad":
        adata = sc.read_h5ad(path)
        logger.info(f"Loaded H5AD: {adata.n_obs:,} × {adata.n_vars:,}")

    elif path.suffix in (".h5", ".hdf5"):
        adata = sc.read_10x_h5(path)
        logger.info(f"Loaded 10x H5: {adata.n_obs:,} × {adata.n_vars:,}")

    elif path.is_dir():
        for subdir in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix", ""):
            candidate = path / subdir if subdir else path
            if (candidate / "matrix.mtx.gz").exists() or (candidate / "matrix.mtx").exists():
                adata = sc.read_10x_mtx(candidate, var_names=var_names, cache=True)
                logger.info(f"Loaded 10x MTX: {adata.n_obs:,} × {adata.n_vars:,}")
                break
        else:
            raise FileNotFoundError(
                f"No recognisable 10x matrix found under {path}. "
                "Expected matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz."
            )

    else:
        raise ValueError(
            f"Unsupported format: {path.suffix}. "
            "Provide a 10x directory, .h5, or .h5ad file."
        )

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    if sample_name is not None:
        adata.obs["sample"] = sample_name

    return adata


def save_adata(adata: sc.AnnData, path: Union[str, Path]) -> None:
    """Save AnnData to compressed H5AD."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression="gzip")
    logger.info(f"Saved {adata.n_obs:,} × {adata.n_vars:,} → {path}")


def pretty_counts(adata: sc.AnnData) -> str:
    """Return a short human-readable summary of cell and gene counts."""
    return f"{adata.n_obs:,} cells × {adata.n_vars:,} genes"


# ── QC ────────────────────────────────────────────────────────────────────────

def add_qc_metrics(adata: sc.AnnData, species: str = "human") -> sc.AnnData:
    """
    Annotate genes and calculate per-cell QC metrics in-place.

    Adds to adata.var:
        mt, ribo, hb  (boolean gene class flags)

    Adds to adata.obs (via scanpy calculate_qc_metrics):
        total_counts, n_genes_by_counts, pct_counts_mt, pct_counts_ribo, pct_counts_hb
    """
    if species.lower() in ("human", "hs", "homo_sapiens"):
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.contains(r"^HB[^PS]", regex=True)
    elif species.lower() in ("mouse", "mm", "mus_musculus"):
        adata.var["mt"] = adata.var_names.str.startswith("mt-")
        adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
        adata.var["hb"] = adata.var_names.str.contains(r"^Hb[^ps]", regex=True)
    else:
        logger.warning(f"Unknown species '{species}'. Skipping gene-class annotation.")
        adata.var["mt"] = False
        adata.var["ribo"] = False
        adata.var["hb"] = False

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[20, 50, 100, 200, 500],
        log1p=True,
        inplace=True,
    )
    return adata


def mad_outlier(
    values: Union[np.ndarray, pd.Series],
    n_mads: float = 5.0,
    log_transform: bool = False,
) -> np.ndarray:
    """
    Return boolean mask of outliers detected by Median Absolute Deviation.

    Parameters
    ----------
    values : array-like
        Numeric values to test.
    n_mads : float
        Number of MADs around the median to define outlier boundary.
    log_transform : bool
        Log1p-transform before computing MAD (useful for count data).

    Returns
    -------
    Boolean array; True where value is an outlier.
    """
    v = np.asarray(values, dtype=float)
    if log_transform:
        v = np.log1p(v)
    med = np.median(v)
    mad = np.median(np.abs(v - med))
    if mad == 0:
        mad = np.mean(np.abs(v - med))  # fallback
    lower = med - n_mads * mad
    upper = med + n_mads * mad
    return (v < lower) | (v > upper)


def adaptive_qc_filter(
    adata: sc.AnnData,
    n_mads: float = 5.0,
    max_mt_pct: Optional[float] = None,
) -> pd.Series:
    """
    Compute per-cell QC pass/fail using MAD-based adaptive thresholds.

    Flags as outliers: extreme total counts, extreme gene counts, high MT%.
    An absolute cap on MT% can be layered on top (recommended: 20-25%).

    Returns
    -------
    pd.Series[bool] — True means the cell passes QC.
    """
    obs = adata.obs
    out_counts = mad_outlier(obs["total_counts"], n_mads, log_transform=True)
    out_genes = mad_outlier(obs["n_genes_by_counts"], n_mads, log_transform=True)
    out_mt = mad_outlier(obs["pct_counts_mt"], n_mads, log_transform=False)

    pass_mask = ~(out_counts | out_genes | out_mt)
    if max_mt_pct is not None:
        pass_mask &= obs["pct_counts_mt"] <= max_mt_pct

    return pd.Series(pass_mask.values, index=obs.index)


# ── Clustering ────────────────────────────────────────────────────────────────

def leiden_sweep(
    adata: sc.AnnData,
    resolutions: List[float],
    key_prefix: str = "leiden",
    random_state: int = 42,
) -> sc.AnnData:
    """
    Run Leiden clustering at multiple resolutions and store each in adata.obs.

    Adds columns {key_prefix}_{res} for each resolution.
    """
    for res in resolutions:
        key = f"{key_prefix}_{res:.2f}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=random_state)
        logger.info(f"Resolution {res}: {adata.obs[key].nunique()} clusters")
    return adata


# ── Pseudobulk Aggregation ────────────────────────────────────────────────────

def pseudobulk(
    adata: sc.AnnData,
    sample_col: str = "sample",
    group_col: str = "cell_type",
    layer: Optional[str] = None,
    min_cells: int = 10,
) -> pd.DataFrame:
    """
    Aggregate single-cell counts into pseudobulk samples for each group.

    Sums raw counts per (sample × cell_type) combination.
    Only combinations with >= min_cells contribute.

    Returns
    -------
    pd.DataFrame with shape (genes × pseudobulk_samples).
    Columns are named "{sample}__{cell_type}".
    """
    X = adata.layers[layer] if layer else adata.X
    if issparse(X):
        X = X.toarray()

    results = {}
    for grp, idx in adata.obs.groupby([sample_col, group_col]).groups.items():
        if len(idx) < min_cells:
            continue
        col_name = f"{grp[0]}__{grp[1]}"
        results[col_name] = X[adata.obs.index.get_indexer(idx), :].sum(axis=0)

    df = pd.DataFrame(results, index=adata.var_names)
    return df


# ── Marker Genes ──────────────────────────────────────────────────────────────

def get_marker_genes(organism: str = "human") -> Dict[str, List[str]]:
    """
    Canonical marker gene sets for common cell types.

    Covers immune (T, B, NK, Myeloid, DC, Plasma), stromal
    (Fibroblast, Endothelial), and epithelial populations.

    Parameters
    ----------
    organism : str
        "human" or "mouse".
    """
    if organism.lower() in ("human", "hs", "homo_sapiens"):
        return {
            "T cell": ["CD3D", "CD3E", "CD3G", "TRAC"],
            "CD4+ T cell": ["CD3D", "CD4", "IL7R", "TCF7"],
            "CD8+ T cell": ["CD3D", "CD8A", "CD8B", "GZMK"],
            "Regulatory T cell": ["CD3D", "CD4", "FOXP3", "IL2RA", "CTLA4"],
            "Exhausted T cell": ["CD3D", "PDCD1", "LAG3", "HAVCR2", "TIGIT"],
            "NK cell": ["NCAM1", "NKG7", "GNLY", "KLRD1", "KLRB1"],
            "B cell": ["CD19", "MS4A1", "CD79A", "CD79B", "BANK1"],
            "Plasma cell": ["MZB1", "IGHG1", "CD38", "SDC1", "XBP1"],
            "Classical monocyte": ["CD14", "LYZ", "CST3", "FCN1", "S100A8"],
            "Non-classical monocyte": ["FCGR3A", "CDKN1C", "VMO1", "CX3CR1"],
            "cDC1": ["CLEC9A", "XCR1", "CADM1"],
            "cDC2": ["FCER1A", "CLEC10A", "CD1C", "FCGR2B"],
            "pDC": ["CLEC4C", "IL3RA", "LILRA4", "JCHAIN"],
            "Macrophage": ["CD68", "CD163", "MRC1", "APOE", "C1QA"],
            "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9", "CXCR2"],
            "Mast cell": ["KIT", "CPA3", "TPSAB1", "TPSB2"],
            "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
            "Fibroblast": ["DCN", "COL1A1", "COL1A2", "LUM", "FAP"],
            "Endothelial": ["PECAM1", "VWF", "CDH5", "SELE", "ENG"],
            "Pericyte": ["ACTA2", "RGS5", "MYH11", "PDGFRB"],
            "Smooth muscle": ["ACTA2", "MYH11", "CNN1", "TAGLN"],
        }
    elif organism.lower() in ("mouse", "mm", "mus_musculus"):
        return {
            "T cell": ["Cd3d", "Cd3e", "Cd3g", "Trac"],
            "CD4+ T cell": ["Cd3d", "Cd4", "Il7r", "Tcf7"],
            "CD8+ T cell": ["Cd3d", "Cd8a", "Cd8b1", "Gzmk"],
            "Regulatory T cell": ["Cd3d", "Cd4", "Foxp3", "Il2ra"],
            "NK cell": ["Ncam1", "Nkg7", "Gzma", "Klrd1"],
            "B cell": ["Cd19", "Ms4a1", "Cd79a", "Cd79b"],
            "Plasma cell": ["Mzb1", "Cd38", "Sdc1", "Xbp1"],
            "Monocyte": ["Cd14", "Lyz2", "Csf1r", "Fcgr3"],
            "Macrophage": ["Cd68", "Cd163", "Mrc1", "Apoe", "C1qa"],
            "Dendritic cell": ["H2-Ab1", "Itgax", "Siglech"],
            "Neutrophil": ["S100a8", "S100a9", "Csf3r", "Ly6g"],
            "Fibroblast": ["Dcn", "Col1a1", "Col1a2", "Pdgfra"],
            "Endothelial": ["Pecam1", "Cdh5", "Tie1", "Vwf"],
            "Epithelial": ["Epcam", "Krt8", "Krt18", "Cdh1"],
            "Pericyte": ["Acta2", "Rgs5", "Pdgfrb"],
        }
    else:
        raise ValueError(f"Unknown organism '{organism}'. Use 'human' or 'mouse'.")
