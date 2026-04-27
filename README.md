# sc-workflow

**A complete, reproducible single-cell RNA-seq analysis pipeline — from raw FASTQ files to advanced multi-modal analyses.**

Available in both **Python** ([scanpy](https://scanpy.readthedocs.io)) and **R** ([Seurat](https://seuratproject.org) / Bioconductor). Each step saves an intermediate file (`.h5ad` for Python, `.rds` for R) that feeds directly into the next, giving you a fully auditable, reproducible chain of evidence from raw counts to biological insight.

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

| Script | Language | Purpose | Input | Output |
|--------|----------|---------|-------|--------|
| `bash/00_fastq_to_counts.sh` | Bash | FASTQ download (SRA), FastQC/MultiQC, Cell Ranger | FASTQ / SRA accession | `data/raw/<sample>/outs/` |
| `Python/01_qc_filtering.py` | Python | QC metrics, doublets, adaptive cell filtering | Cell Ranger MTX / `.h5ad` | `results/01_qc/adata_qc.h5ad` |
| `Python/02_integration_clustering.py` | Python | Normalise, HVG, PCA, Harmony, UMAP, Leiden | `adata_qc.h5ad` | `results/02_integration/adata_integrated.h5ad` |
| `Python/03_cell_annotation.py` | Python | CellTypist, marker scoring, manual overrides | `adata_integrated.h5ad` | `results/03_annotation/adata_annotated.h5ad` |
| `Python/04_differential_expression.py` | Python | Pseudobulk DESeq2, proportion testing, Enrichr | `adata_annotated.h5ad` | `results/04_de/` |
| `Python/05_trajectory_analysis.py` | Python | CellRank + diffusion pseudotime / scVelo | `adata_annotated.h5ad` | `results/05_trajectory/` |
| `Python/06_cell_communication.py` | Python | LIANA aggregate (CellChat, NATMI, logFC…) | `adata_annotated.h5ad` | `results/06_communication/` |
| `Python/07_cnv_analysis.py` | Python | inferCNVpy, tumour vs normal classification | `adata_annotated.h5ad` | `results/07_cnv/adata_cnv.h5ad` |
| `Python/08_coexpression_networks.py` | Python | PyWGCNA metacell co-expression modules | `adata_annotated.h5ad` | `results/08_coexpression/` |

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

`bash/00_fastq_to_counts.sh` requires FastQC, MultiQC, SRA Toolkit, and samtools. These are available as compiled binaries via bioconda and are installed into a separate pixi environment called `genomics`.

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
pixi run python Python/01_qc_filtering.py

# Option B — enter the pixi shell (PATH is set automatically)
pixi shell
python Python/01_qc_filtering.py
exit
```

### Using the built-in task shortcuts

```bash
pixi run qc            # Python/01_qc_filtering.py
pixi run integrate     # Python/02_integration_clustering.py
pixi run annotate      # Python/03_cell_annotation.py
pixi run de            # Python/04_differential_expression.py
pixi run trajectory    # Python/05_trajectory_analysis.py
pixi run communication # Python/06_cell_communication.py
pixi run cnv           # Python/07_cnv_analysis.py
pixi run coexpression  # Python/08_coexpression_networks.py
```

### Running the full pipeline in order

```bash
pixi run pipeline
```

> Each task depends on the previous one completing successfully; the pipeline stops on any error.

---

## Detailed Usage

### Script 00 — FASTQ to Count Matrix *(bash, genomics environment)*

```bash
# Edit the CONFIGURATION block at the top of bash/00_fastq_to_counts.sh, then:
chmod +x bash/00_fastq_to_counts.sh
pixi shell -e genomics
bash bash/00_fastq_to_counts.sh
```

Key configuration variables:

```bash
SPECIES="human"                          # or "mouse"
SAMPLE_IDS=(      "pbmc_ctrl"  )
SRA_ACCESSIONS=(  "SRR12345678")         # empty string "" if FASTQs already on disk
SAMPLE_PREFIXES=( "pbmc_ctrl"  )
N_CORES=16
MEM_GB=128
```

### Scripts 01–08 — Python Analysis Pipeline

Each script follows the same pattern:

1. Open the script in any editor.
2. Edit the `CONFIGURATION` block near the top.
3. Run via `pixi run <task>` or `pixi run python Python/NN_*.py`.

**`Python/01_qc_filtering.py` — QC & Filtering**

```python
DATA_PATH   = "data/raw/sample1/outs/filtered_feature_bc_matrix"
SAMPLE_NAME = "sample1"
SPECIES     = "human"    # or "mouse"
MIN_GENES   = 500
MAX_MT_PCT  = 20
OUTPUT_PATH = "results/01_qc/adata_qc.h5ad"
```

**`Python/02_integration_clustering.py` — Integration & Clustering**

```python
INPUT_PATHS            = ["results/01_qc/adata_qc.h5ad"]
BATCH_KEY              = "sample"
N_HVGS                 = 3000
CLUSTERING_RESOLUTIONS = [0.3, 0.5, 0.8, 1.0]
OUTPUT_PATH            = "results/02_integration/adata_integrated.h5ad"
```

**`Python/03_cell_annotation.py` — Cell Annotation**

```python
INPUT_PATH         = "results/02_integration/adata_integrated.h5ad"
CELLTYPIST_MODEL   = "Immune_All_High.pkl"   # auto-downloads if missing
ORGANISM           = "human"
MANUAL_ANNOTATION  = {}   # e.g. {"5": "Plasmablast"}
OUTPUT_PATH        = "results/03_annotation/adata_annotated.h5ad"
```

**`Python/04_differential_expression.py` — Differential Expression**

```python
REFERENCE_CONDITION = "Control"
TEST_CONDITION      = "Treatment"
PADJ_CUTOFF         = 0.05
LFC_CUTOFF          = 0.5
OUTPUT_DIR          = "results/04_de"
```

**Scripts 05–08** follow the same pattern. Refer to the `CONFIGURATION` block at the top of each script for the full parameter list.

---

## R Workflow

An equivalent pipeline is available in `R/` using idiomatic R/Bioconductor packages. It produces the same results as the Python workflow — choose whichever language you prefer.

| Script | Language | R packages | Python equivalent |
|--------|----------|-----------|-------------------|
| `bash/00_fastq_to_counts.sh` | Bash | — shared with Python workflow | — shared script |
| `R/01_qc_filtering.R` | R | Seurat, scDblFinder | scanpy, scrublet |
| `R/02_integration_clustering.R` | R | Seurat, harmony | scanpy, harmonypy |
| `R/03_cell_annotation.R` | R | SingleR, celldex | celltypist |
| `R/04_differential_expression.R` | R | DESeq2, enrichR | pydeseq2, gseapy |
| `R/05_trajectory_analysis.R` | R | monocle3 | cellrank |
| `R/06_cell_communication.R` | R | CellChat | liana-py |
| `R/07_cnv_analysis.R` | R | infercnv | infercnvpy |
| `R/08_coexpression_networks.R` | R | WGCNA | PyWGCNA |

Data flow uses `.rds` (Seurat objects) instead of `.h5ad`. The quantification step (`bash/00_fastq_to_counts.sh`) is a **shared bash script** — run it once and its count matrix output feeds both the Python and R analysis pipelines. Scripts 01–08 are pure R: edit the `CONFIGURATION` block at the top, then run with `Rscript`.

---

### R Requirements

| Requirement | Details |
|-------------|---------|
| **R** | ≥ 4.4.0 |
| **OS** | Linux (recommended), macOS, Windows (scripts 01–08) |
| **RAM** | ≥ 32 GB recommended |
| **Disk** | ≥ 50 GB for `.rds` intermediates |

---

### R — Step 1: Install R

Download and install R ≥ 4.4.0 from [CRAN](https://cran.r-project.org/).

**Linux system libraries** (required before installing R packages):

```bash
# Ubuntu / Debian
sudo apt-get install -y \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libhdf5-dev libfontconfig1-dev libfreetype6-dev \
  libharfbuzz-dev libfribidi-dev libpng-dev libtiff-dev

# Fedora / RHEL
sudo dnf install -y \
  libcurl-devel openssl-devel libxml2-devel \
  hdf5-devel fontconfig-devel freetype-devel
```

**macOS** (Homebrew):

```bash
brew install hdf5 libgit2 pkg-config
```

---

### R — Step 2: Install pak

[pak](https://pak.r-lib.org) is a fast R package manager that installs from CRAN, Bioconductor, and GitHub with a single unified function — no separate `BiocManager::install()` calls needed.

Run this **once** in your R console:

```r
install.packages("pak")
```

Verify:

```r
packageVersion("pak")
```

---

### R — Step 3: Install all R packages

Run the following block in your R console (or save as `install_packages.R` and run `Rscript install_packages.R`). pak resolves all dependencies, downloads in parallel, and handles CRAN, Bioconductor, and GitHub automatically.

```r
pak::pak(c(
  # ── Core ─────────────────────────────────────────────────────────
  "Seurat",
  "SeuratObject",
  "ggplot2",
  "patchwork",
  "dplyr",
  "ggrepel",
  "Matrix",
  "FNN",
  "tools",

  # ── QC & Doublet Detection ────────────────────────────────────────
  "bioc::scDblFinder",
  "bioc::SingleCellExperiment",
  "bioc::BiocParallel",

  # ── Batch Correction ──────────────────────────────────────────────
  "harmony",

  # ── Cell Type Annotation ──────────────────────────────────────────
  "bioc::SingleR",
  "bioc::celldex",

  # ── Differential Expression ───────────────────────────────────────
  "bioc::DESeq2",
  "enrichR",

  # ── Trajectory Analysis ───────────────────────────────────────────
  "cole-trapnell-lab/monocle3",

  # ── Cell-Cell Communication ───────────────────────────────────────
  "jinworks/CellChat",

  # ── CNV Analysis ──────────────────────────────────────────────────
  "bioc::infercnv",

  # ── Co-expression Networks ────────────────────────────────────────
  "WGCNA"
))
```

> **Note:** `monocle3` and `CellChat` are installed from GitHub. pak handles authentication-free public repos automatically.

---

### R — Step 4: Download the inferCNV gene order file *(script 07 only)*

Script 07 requires a chromosome-position file for human or mouse genes:

```bash
# Human (hg38)
wget -P data/raw/ \
  https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt

# Mouse (mm10)
wget -P data/raw/ \
  https://data.broadinstitute.org/Trinity/CTAT/cnv/mm10_gencode_vM25.txt
```

Set the path in `07_cnv_analysis.R`:

```r
GENE_ORDER_FILE <- "data/raw/hg38_gencode_v27.txt"
```

---

### R — Quick Start

**Step 0 — Quantification (shared bash script)**

Edit the `CONFIGURATION` block at the top of `bash/00_fastq_to_counts.sh`, then run inside the pixi genomics environment:

```bash
chmod +x bash/00_fastq_to_counts.sh
pixi shell -e genomics
bash bash/00_fastq_to_counts.sh
```

This downloads FASTQs, runs FastQC/MultiQC, and runs Cell Ranger. The output count matrices feed both the Python and R analysis pipelines.

**Steps 01–08 — R analysis scripts**

Each R script has a `CONFIGURATION` block at the top. Edit the paths and parameters, then run:

```bash
# Run a single script
Rscript R/01_qc_filtering.R

# Run the full analysis pipeline in order
Rscript R/01_qc_filtering.R            && \
Rscript R/02_integration_clustering.R  && \
Rscript R/03_cell_annotation.R         && \
Rscript R/04_differential_expression.R && \
Rscript R/05_trajectory_analysis.R     && \
Rscript R/06_cell_communication.R      && \
Rscript R/07_cnv_analysis.R            && \
Rscript R/08_coexpression_networks.R
```

Or source interactively from an R session / RStudio:

```r
source("R/01_qc_filtering.R")
```

---

### R Script Configuration — Key Parameters

**Script 01 — QC & Filtering**

```r
DATA_PATH   <- "data/raw/sample1/outs/filtered_feature_bc_matrix"
SAMPLE_NAME <- "sample1"
SPECIES     <- "human"      # or "mouse"
MIN_GENES   <- 500
MAX_GENES   <- 8000
MAX_MT_PCT  <- 20
```

**Script 02 — Integration & Clustering**

```r
INPUT_PATHS            <- c("results/01_qc/sample1_qc.rds")
BATCH_KEY              <- "sample"
N_HVGS                 <- 3000L
CLUSTERING_RESOLUTIONS <- c(0.3, 0.5, 0.8, 1.0)
DEFAULT_RESOLUTION     <- 0.5
```

**Script 03 — Cell Annotation**

```r
CLUSTER_KEY        <- "leiden_0.50"
SINGLER_REFERENCE  <- "HumanPrimaryCellAtlasData"  # or "BlueprintEncodeData"
MANUAL_ANNOTATION  <- list()   # e.g. list("0" = "CD4+ T cell", "3" = "Monocyte")
```

**Script 04 — Differential Expression**

```r
CONDITION_COL       <- "condition"
REFERENCE_CONDITION <- "Control"
TEST_CONDITION      <- "Treatment"
PADJ_CUTOFF         <- 0.05
LFC_CUTOFF          <- 0.5
```

**Scripts 05–08** follow the same pattern. Refer to the `CONFIGURATION` block inside each script for the full parameter list.

---

### R Troubleshooting

**`Error: package 'XYZ' is not available`**

Run `pak::pak("XYZ")` individually to see a detailed error. Most failures are missing system libraries — install them with apt/brew (see Step 1 above), then retry.

**pak fails on Bioconductor packages**

pak auto-detects the correct Bioconductor version for your R installation. If it fails, set the version explicitly:

```r
options(pkg.bioc_version = "3.20")   # match your R 4.4.x install
pak::pak("bioc::DESeq2")
```

**CellChat or monocle3 GitHub install times out**

Install them separately with a longer timeout:

```r
options(timeout = 600)
pak::pak("jinworks/CellChat")
pak::pak("cole-trapnell-lab/monocle3")
```

**`Error in JoinLayers(seu)`: Seurat v5 layer issue**

All R scripts target Seurat v5. If you have Seurat v4 installed, update it:

```r
pak::pak("Seurat")
```

**High memory use during WGCNA (script 08)**

Reduce metacell count or restrict to fewer genes:

```r
METACELL_K        <- 10L    # fewer neighbours per metacell
MIN_EXPR_FRACTION <- 0.10   # stricter gene filter (default 0.05)
```

**Script 07 fails: `gene_order_file not found`**

Download the file for your genome version (see Step 4 above) and update `GENE_ORDER_FILE` in the script.

**`00_fastq_to_counts.sh`: `permission denied`**

Make the script executable first:

```bash
chmod +x R/00_fastq_to_counts.sh
```

**`00_fastq_to_counts.sh`: tool not found (fastqc, fasterq-dump, etc.)**

These tools are provided by the pixi `genomics` environment. Activate it before running the bash script:

```bash
pixi shell -e genomics
bash R/00_fastq_to_counts.sh
```

---

## Directory Structure

```
sc-workflow/
├── pixi.toml                   # pixi environment definition (Python)
├── environment.yml             # conda-compatible environment (pixi can import)
├── Python/                     # 9 plain-Python analysis scripts (run in order)
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
├── R/                          # R/bash analysis scripts (run in order)
│   ├── 00_fastq_to_counts.sh   # bash — shared quantification script (same as Python 00)
│   ├── 01_qc_filtering.R
│   ├── 02_integration_clustering.R
│   ├── 03_cell_annotation.R
│   ├── 04_differential_expression.R
│   ├── 05_trajectory_analysis.R
│   ├── 06_cell_communication.R
│   ├── 07_cnv_analysis.R
│   └── 08_coexpression_networks.R
│
├── src/sc_workflow/
│   ├── __init__.py
│   └── utils.py                # Shared Python helpers
│
├── data/
│   └── raw/                    # Cell Ranger output / downloaded FASTQs
│       └── fastqs/
│
└── results/                    # Auto-created by both Python and R scripts
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

| Analysis step | R package | Python package |
|---------------|-----------|----------------|
| Single-cell data structure | Seurat / SeuratObject | scanpy + anndata |
| Doublet detection | scDblFinder | scrublet / doubletdetection |
| Batch correction | harmony | harmonypy |
| Leiden clustering | `FindClusters(algorithm=4)` | `sc.tl.leiden` |
| Cell type annotation | SingleR + celldex | celltypist |
| Differential expression | DESeq2 | pydeseq2 |
| Trajectory analysis | monocle3 | cellrank |
| Cell communication | CellChat | liana-py |
| CNV inference | infercnv | infercnvpy |
| Co-expression networks | WGCNA | PyWGCNA |
| Pathway enrichment | enrichR | gseapy (Enrichr) + decoupler-py |

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
