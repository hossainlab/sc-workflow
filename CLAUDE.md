# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A Python-based single-cell RNA-seq analysis workflow implemented as interactive **Marimo** notebooks.  Covers the complete pipeline from raw counts to advanced analyses, designed for laboratory research use.

## Commands

```bash
# Install all dependencies
pip install -r requirements.txt

# Run a notebook interactively (edit mode)
marimo edit notebooks/01_qc_filtering.py

# Run a notebook as an app (no code visible)
marimo run notebooks/01_qc_filtering.py

# Install the sc_workflow utility package in editable mode
pip install -e src/
```

## Analysis Scripts (run in order)

All scripts are **plain Python** — no Marimo. Edit the `CONFIGURATION` section at
the top of each file, then run with `python notebooks/<script>.py`.

| Script | Purpose | Input | Output |
|---|---|---|---|
| `00_fastq_to_counts.py` | FASTQ download (SRA), FastQC/MultiQC, Cell Ranger count + aggr | Raw FASTQ / SRA accession | `data/raw/<sample>/outs/filtered_feature_bc_matrix/` |
| `01_qc_filtering.py` | QC metrics, Scrublet doublets, cell/gene filtering | Cell Ranger output / .h5ad | `results/01_qc/adata_qc.h5ad` |
| `02_integration_clustering.py` | Normalise, HVG, PCA, Harmony, UMAP, Leiden | `adata_qc.h5ad` | `results/02_integration/adata_integrated.h5ad` |
| `03_cell_annotation.py` | CellTypist, marker scoring, manual override dict | `adata_integrated.h5ad` | `results/03_annotation/adata_annotated.h5ad` |
| `04_differential_expression.py` | Pseudobulk DESeq2, proportion testing, Enrichr | `adata_annotated.h5ad` | `results/04_de/` |
| `05_trajectory_analysis.py` | CellRank + diffusion pseudotime / scVelo | `adata_annotated.h5ad` | `results/05_trajectory/` |
| `06_cell_communication.py` | LIANA aggregate (CellChat, NATMI, logFC…) | `adata_annotated.h5ad` | `results/06_communication/` |
| `07_cnv_analysis.py` | inferCNVpy, tumour vs normal classification | `adata_annotated.h5ad` | `results/07_cnv/adata_cnv.h5ad` |
| `08_coexpression_networks.py` | PyWGCNA metacell co-expression modules | `adata_annotated.h5ad` | `results/08_coexpression/` |

## Architecture

```
sc-workflow/
├── requirements.txt            # All Python (pip) dependencies
├── src/sc_workflow/
│   ├── __init__.py
│   └── utils.py                # Shared helpers (load_data, add_qc_metrics,
│                               #  mad_outlier, pseudobulk, leiden_sweep, etc.)
├── notebooks/                  # 9 Marimo notebooks (00–08)
│   └── 00_fastq_to_counts.py  # FASTQ → count matrix (bash/Cell Ranger)
├── data/raw/                   # Cell Ranger output lands here
│   └── fastqs/                 # FASTQ files (SRA downloads or local)
└── results/                    # Auto-created by notebooks 01–08
```

## System Dependencies (not in requirements.txt)

Notebook 00 requires system-level tools installed separately:

| Tool | Install |
|---|---|
| Cell Ranger ≥ 8.0 | https://www.10xgenomics.com/support/software/cell-ranger/downloads |
| FastQC | `conda install -c bioconda fastqc` |
| MultiQC | `conda install -c bioconda multiqc` |
| SRA Toolkit (fasterq-dump) | `conda install -c bioconda sra-tools` |
| samtools | `conda install -c bioconda samtools` |

## Key Design Decisions

**Marimo reactive model:** Each notebook cell is a pure function. Variables are shared via return values; UI elements (`mo.ui.*`) trigger automatic re-computation of dependent cells when changed.  Use `mo.stop(condition, output)` to halt a cell on unmet conditions.

**Data flow:** Notebooks communicate through `.h5ad` files.  Every notebook loads the previous step's output, adds annotations/results to `adata.obs` / `adata.uns` / `adata.obsm`, and saves a new `.h5ad`.

**Raw counts:** Always preserved in `adata.layers["counts"]` (added in notebook 02).  Required by DESeq2 pseudobulk (notebook 04) and inferCNVpy (notebook 07).

**Per-sample processing:** QC (notebook 01) must be run once per sample.  Normalisation and HVG selection in notebook 02 are performed per sample before merging.

**Pseudobulk for DE:** The `pseudobulk()` utility in `utils.py` aggregates raw counts per (sample × cell_type). Never use direct cell-level Wilcoxon results as final published DE findings.

**Species support:** `add_qc_metrics()` and `get_marker_genes()` both accept `species="human"` or `species="mouse"`.

## Shared Utilities (`src/sc_workflow/utils.py`)

| Function | Purpose |
|---|---|
| `load_data(path)` | Loads 10x MTX dir, .h5, or .h5ad |
| `add_qc_metrics(adata, species)` | Annotates mt/ribo/hb genes, calls `sc.pp.calculate_qc_metrics` |
| `mad_outlier(values, n_mads)` | MAD-based outlier detection for adaptive QC thresholds |
| `adaptive_qc_filter(adata, n_mads)` | Returns per-cell pass/fail Series |
| `pseudobulk(adata, sample_col, group_col)` | Sums raw counts per sample×celltype |
| `leiden_sweep(adata, resolutions)` | Runs Leiden at multiple resolutions |
| `get_marker_genes(organism)` | Returns canonical marker gene dict |

## Python ↔ R Tool Mapping

| R tool (tutorials) | Python equivalent used here |
|---|---|
| Seurat | scanpy + anndata |
| scDblFinder | scrublet |
| Harmony | harmonypy |
| Leiden algorithm | `sc.tl.leiden` (leidenalg) |
| SingleR / scType | celltypist |
| DESeq2 | pydeseq2 |
| Monocle 3 / Slingshot | cellrank |
| CellChat / NicheNet | liana-py |
| CopyKAT / inferCNV | infercnvpy |
| hdWGCNA | PyWGCNA |
