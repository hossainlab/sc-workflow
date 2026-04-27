# sc-workflow

**A complete, reproducible single-cell RNA-seq analysis pipeline — from raw FASTQ files to advanced multi-modal analyses.**

Built with Python, [scanpy](https://scanpy.readthedocs.io), and [Marimo](https://marimo.io) reactive notebooks. Each step saves an annotated `.h5ad` file that feeds directly into the next, giving you a fully auditable, reproducible chain of evidence from raw counts to biological insight.

---

## Pipeline at a Glance

```
FASTQ / SRA
    │
    ▼
00  FASTQ → Counts      Cell Ranger count + aggr, FastQC, MultiQC
    │
    ▼
01  QC & Filtering       Scrublet doublets, MAD-based adaptive thresholds
    │
    ▼
02  Integration          Normalise, HVG, PCA, Harmony, UMAP, Leiden
    │
    ▼
03  Cell Annotation      CellTypist, canonical markers, manual overrides
    │
    ├──▶ 04  Differential Expression   Pseudobulk DESeq2, Enrichr pathway
    ├──▶ 05  Trajectory Analysis       CellRank + scVelo RNA velocity
    ├──▶ 06  Cell Communication        LIANA aggregate (CellChat, NATMI…)
    ├──▶ 07  CNV Analysis              inferCNVpy tumour vs normal
    └──▶ 08  Co-expression Networks    PyWGCNA metacell modules
```

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `00_fastq_to_counts.py` | FASTQ download (SRA), FastQC/MultiQC, Cell Ranger | FASTQ / SRA accession | `data/raw/<sample>/outs/` |
| `01_qc_filtering.py` | QC metrics, doublets, adaptive cell filtering | Cell Ranger MTX / `.h5ad` | `results/01_qc/adata_qc.h5ad` |
| `02_integration_clustering.py` | Normalise, HVG, PCA, Harmony, UMAP, Leiden | `adata_qc.h5ad` | `results/02_integration/adata_integrated.h5ad` |
| `03_cell_annotation.py` | CellTypist, marker scoring, manual overrides | `adata_integrated.h5ad` | `results/03_annotation/adata_annotated.h5ad` |
| `04_differential_expression.py` | Pseudobulk DESeq2, proportion testing, Enrichr | `adata_annotated.h5ad` | `results/04_de/` |
| `05_trajectory_analysis.py` | CellRank + diffusion pseudotime / scVelo | `adata_annotated.h5ad` | `results/05_trajectory/` |
| `06_cell_communication.py` | LIANA aggregate (CellChat, NATMI, logFC…) | `adata_annotated.h5ad` | `results/06_communication/` |
| `07_cnv_analysis.py` | inferCNVpy, tumour vs normal classification | `adata_annotated.h5ad` | `results/07_cnv/adata_cnv.h5ad` |
| `08_coexpression_networks.py` | PyWGCNA metacell co-expression modules | `adata_annotated.h5ad` | `results/08_coexpression/` |

---

## Requirements

| Requirement | Details |
|-------------|---------|
| **OS** | Linux (recommended), macOS (Intel or Apple Silicon), Windows (scripts 01–08 only) |
| **Python** | 3.11 or 3.12 |
| **RAM** | ≥ 32 GB recommended (16 GB minimum for small datasets) |
| **Disk** | ≥ 50 GB free for intermediate `.h5ad` files + Cell Ranger output |
| **GPU** | Optional — used by scVI-tools and scVelo if a CUDA device is present |

> **Cell Ranger** (script 00) requires Linux x86-64 and a 10x Genomics licence. Scripts 01–08 run on all platforms.

---

## Installation

### 1 — Install pixi

[pixi](https://pixi.sh) is a fast, cross-platform package manager that handles both conda packages and PyPI packages in a single lock file. No conda required.

**Linux / macOS**

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

**Windows (PowerShell)**

```powershell
iwr -useb https://pixi.sh/install.ps1 | iex
```

After installation, open a new terminal (or run `source ~/.bashrc` / restart PowerShell) so `pixi` is on your `PATH`.

Verify:

```bash
pixi --version
```

---

### 2 — Clone the repository

```bash
git clone https://github.com/JHossain/sc-workflow.git
cd sc-workflow
```

---

### 3 — Install the default environment

This installs all Python analysis packages (scripts 01–08) into an isolated environment managed by pixi.

```bash
pixi install
```

The first run resolves and downloads all packages. Subsequent runs use the cached lock file (`pixi.lock`) and are nearly instant.

---

### 4 — Install the genomics environment *(optional — Linux/macOS only)*

Script `00_fastq_to_counts.py` requires FastQC, MultiQC, SRA Toolkit, and samtools. These are available as compiled binaries via bioconda and are installed into a separate pixi environment called `genomics`.

```bash
pixi install -e genomics
```

> **Windows users:** Skip this step. System genomics tools are not available on Windows. You can still run scripts 01–08 on count matrices you already have.

---

### 5 — Install Cell Ranger *(script 00 only)*

Cell Ranger must be installed separately because it requires a free 10x Genomics registration.

1. Download from [10x Genomics Software Downloads](https://www.10xgenomics.com/support/software/cell-ranger/downloads) (version ≥ 8.0).
2. Extract and add to `PATH`:

```bash
export PATH=/path/to/cellranger-8.x.x:$PATH
```

3. Verify:

```bash
cellranger --version
```

---

### Alternative: import environment.yml

If you prefer the conda environment file format (e.g. to import into an existing workflow), an `environment.yml` is also provided. You can import it into pixi:

```bash
pixi init --import environment.yml
pixi install
```

---

## Quick Start

### Running a single script

Each script has a `CONFIGURATION` block at the top. Open the script, edit the paths and parameters, then run it inside the pixi environment:

```bash
# Option A — prefix every command with pixi run
pixi run python scripts/01_qc_filtering.py

# Option B — enter the pixi shell (PATH is set automatically)
pixi shell
python scripts/01_qc_filtering.py
exit
```

### Using the built-in task shortcuts

```bash
pixi run qc            # scripts/01_qc_filtering.py
pixi run integrate     # scripts/02_integration_clustering.py
pixi run annotate      # scripts/03_cell_annotation.py
pixi run de            # scripts/04_differential_expression.py
pixi run trajectory    # scripts/05_trajectory_analysis.py
pixi run communication # scripts/06_cell_communication.py
pixi run cnv           # scripts/07_cnv_analysis.py
pixi run coexpression  # scripts/08_coexpression_networks.py
```

### Running the full pipeline in order

```bash
pixi run pipeline
```

> Each task depends on the previous one completing successfully; the pipeline stops on any error.

---

## Detailed Usage

### Script 00 — FASTQ to Count Matrix *(genomics environment)*

```bash
# Edit CONFIGURATION in scripts/00_fastq_to_counts.py, then:
pixi run -e genomics python scripts/00_fastq_to_counts.py
```

Key configuration options:

```python
CONFIGURATION = {
    "sra_accessions": ["SRR12345678"],   # leave empty to skip SRA download
    "fastq_dir":      "data/raw/fastqs",
    "genome_ref":     "/path/to/refdata-gex-GRCh38-2024-A",
    "output_dir":     "data/raw",
    "n_cores":        16,
}
```

### Scripts 01–08 — Analysis Pipeline

Each script follows the same pattern:

1. Open the script in any editor.
2. Edit the `CONFIGURATION` dict near the top.
3. Run via `pixi run <task>` or `pixi run python scripts/NN_*.py`.

**Script 01 — QC & Filtering**

```python
CONFIGURATION = {
    "input_path": "data/raw/sample1/outs/filtered_feature_bc_matrix",
    "species":    "human",   # or "mouse"
    "n_mads":     5,         # MAD multiplier for adaptive thresholds
    "min_genes":  200,
    "output_dir": "results/01_qc",
}
```

**Script 02 — Integration & Clustering**

```python
CONFIGURATION = {
    "input_path":       "results/01_qc/adata_qc.h5ad",
    "batch_key":        "sample",          # column in adata.obs
    "n_hvgs":           3000,
    "leiden_resolution": 0.5,
    "output_dir":       "results/02_integration",
}
```

**Script 03 — Cell Annotation**

```python
CONFIGURATION = {
    "input_path":        "results/02_integration/adata_integrated.h5ad",
    "celltypist_model":  "Immune_All_High.pkl",  # auto-downloads if missing
    "species":           "human",
    "manual_overrides":  {},   # e.g. {"cluster_5": "Plasmablast"}
    "output_dir":        "results/03_annotation",
}
```

**Scripts 04–08** follow the same pattern with their own `CONFIGURATION` blocks. Refer to the comments inside each script for the full list of parameters.

---

## Directory Structure

```
sc-workflow/
├── pixi.toml                   # pixi environment definition (primary)
├── environment.yml             # conda-compatible environment (pixi can import)
├── requirements.txt            # pip-only dependency list (for reference)
│
├── scripts/                    # 9 plain-Python analysis scripts (run in order)
│   ├── 00_fastq_to_counts.py
│   ├── 01_qc_filtering.py
│   ├── 02_integration_clustering.py
│   ├── 03_cell_annotation.py
│   ├── 04_differential_expression.py
│   ├── 05_trajectory_analysis.py
│   ├── 06_cell_communication.py
│   ├── 07_cnv_analysis.py
│   └── 08_coexpression_networks.py
│
├── src/sc_workflow/
│   ├── __init__.py
│   └── utils.py                # Shared helpers (see below)
│
├── data/
│   └── raw/                    # Cell Ranger output / downloaded FASTQs
│       └── fastqs/
│
└── results/                    # Auto-created by scripts 01–08
    ├── 01_qc/
    ├── 02_integration/
    ├── 03_annotation/
    ├── 04_de/
    ├── 05_trajectory/
    ├── 06_communication/
    ├── 07_cnv/
    └── 08_coexpression/
```

---

## Shared Utilities (`src/sc_workflow/utils.py`)

All scripts import from this package, which is installed in editable mode by pixi.

| Function | Signature | Purpose |
|----------|-----------|---------|
| `load_data` | `(path)` | Loads 10x MTX directory, `.h5`, or `.h5ad` |
| `add_qc_metrics` | `(adata, species)` | Annotates MT/ribo/Hb genes, calls `sc.pp.calculate_qc_metrics` |
| `mad_outlier` | `(values, n_mads)` | Returns boolean mask of MAD-based outliers |
| `adaptive_qc_filter` | `(adata, n_mads)` | Returns per-cell pass/fail Series using `mad_outlier` |
| `pseudobulk` | `(adata, sample_col, group_col)` | Sums raw counts per sample × cell-type for DESeq2 |
| `leiden_sweep` | `(adata, resolutions)` | Runs Leiden clustering at multiple resolutions |
| `get_marker_genes` | `(organism)` | Returns canonical marker gene dictionary |

---

## Key Design Decisions

**Raw counts always preserved.** `adata.layers["counts"]` stores the original integer counts added in script 02. Scripts 04 (DESeq2 pseudobulk) and 07 (inferCNVpy) require this layer.

**Pseudobulk DE only.** Script 04 aggregates counts per sample × cell-type before running DESeq2. Direct cell-level Wilcoxon tests are used only for exploratory marker detection, not for published DE results.

**Adaptive QC thresholds.** Script 01 uses MAD-based outlier detection (`n_mads`, default 5) rather than fixed cutoffs. This adapts automatically to dataset-specific distributions.

**Per-sample normalisation.** HVG selection and normalisation in script 02 are performed per sample before merging, avoiding batch effects from joint normalisation.

**Species support.** `add_qc_metrics()` and `get_marker_genes()` both accept `species="human"` or `species="mouse"`.

---

## Python ↔ R Tool Mapping

Coming from a Seurat/R background? Here is the direct equivalent for each tool used in this workflow.

| R tool | Python equivalent |
|--------|-------------------|
| Seurat | scanpy + anndata |
| scDblFinder | scrublet / doubletdetection |
| Harmony | harmonypy |
| Leiden (leidenalg) | `sc.tl.leiden` |
| SingleR / scType | celltypist |
| DESeq2 | pydeseq2 |
| Monocle 3 / Slingshot | cellrank |
| CellChat / NicheNet | liana-py |
| CopyKAT / inferCNV | infercnvpy |
| hdWGCNA | PyWGCNA |
| clusterProfiler / fgsea | decoupler-py + gseapy |

---

## Troubleshooting

### `pixi: command not found`

Restart your terminal after installation, or run `source ~/.bashrc` (bash) / `source ~/.zshrc` (zsh). On Windows, open a new PowerShell window.

### Packages fail to resolve / lock file conflict

Delete the lock file and re-solve:

```bash
rm pixi.lock
pixi install
```

### `ModuleNotFoundError: sc_workflow`

The local package is installed in editable mode. Make sure you ran `pixi install` from the repository root (the directory containing `pixi.toml`).

### Script fails with `KeyError: 'counts'`

You are running a downstream script (04, 07) on an `.h5ad` file produced before script 02 added the `counts` layer. Re-run script 02 from the beginning to regenerate the integrated file.

### Script 00 fails: `cellranger: command not found`

Cell Ranger is not managed by pixi. Follow the [installation steps above](#5--install-cell-ranger-script-00-only) and ensure the `cellranger` binary is on your `PATH`.

### `ERROR: Could not find a version that satisfies the requirement diffxpy`

`diffxpy` is only available for Python ≤ 3.11 on some platforms. If you encounter this, pin Python to 3.11 in `pixi.toml`:

```toml
[dependencies]
python = "3.11.*"
```

Then re-run `pixi install`.

### High memory usage during integration (script 02)

Reduce the number of highly variable genes:

```python
"n_hvgs": 2000   # default is 3000
```

Or subsample cells for initial exploration:

```python
sc.pp.subsample(adata, n_obs=20000)
```

---

## Pixi Cheat Sheet

```bash
pixi install                    # install / sync the default environment
pixi install -e genomics        # install the genomics environment
pixi run <task>                 # run a named task (see pixi.toml [tasks])
pixi shell                      # activate the environment in current shell
pixi shell -e genomics          # activate the genomics environment
pixi list                       # list installed packages
pixi update                     # update packages to latest allowed versions
pixi clean                      # remove cached environments
```

---

## License

MIT — see [LICENSE](LICENSE).

Copyright © 2026 Jubayer Hossain.
