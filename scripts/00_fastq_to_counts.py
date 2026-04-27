#!/usr/bin/env python3
"""
00_fastq_to_counts.py
=====================
FASTQ → Count Matrix preprocessing pipeline for 10x Chromium scRNA-seq data.

WHAT THIS SCRIPT DOES (runs top-to-bottom in this order):
  Step 1 — Check that all required tools are installed
  Step 2 — Download the Cell Ranger reference genome (only needed once, ~10 GB)
  Step 3 — Download FASTQ files from NCBI SRA (skip if files already exist)
  Step 4 — Run FastQC + MultiQC to check raw read quality
  Step 5 — Run Cell Ranger count for each sample
  Step 6 — (Optional) Merge all samples with cellranger aggr
  Step 7 — Print a QC metrics summary table

HOW TO USE:
  1. Edit the CONFIGURATION section below — that is the only section you need to change.
  2. Open a terminal and run:  python 00_fastq_to_counts.py
  3. When it finishes, pass the output paths to notebook 01 (QC & Filtering).

REQUIRED TOOLS (install once before running this script):
  Cell Ranger  →  https://www.10xgenomics.com/support/software/cell-ranger/downloads
  SRA Toolkit  →  conda install -c bioconda sra-tools
  FastQC       →  conda install -c bioconda fastqc
  MultiQC      →  pip install multiqc
  samtools     →  conda install -c bioconda samtools

"""

import subprocess
import shutil
import sys
from pathlib import Path

import pandas as pd


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

# --- Paths --------------------------------------------------------------------

# Where the Cell Ranger pre-built reference genome is (or will be downloaded to).
# You only need to download this once and can reuse it across projects.
REFERENCE_PATH = Path("~/references/cellranger/refdata-gex-GRCh38-2024-A").expanduser()

# Folder where FASTQ files will be stored.
# If downloading from SRA, files are saved here automatically.
# If you already have FASTQ files, put them here (or point to their location).
FASTQ_DIR = Path("data/raw/fastqs")

# Root folder for Cell Ranger output.
# Each sample gets its own subdirectory: data/raw/<sample_id>/outs/
CELLRANGER_OUT_DIR = Path("data/raw")

# --- Reference genome ---------------------------------------------------------

# Which organism are your cells from?
# "human" downloads GRCh38-2024-A  (~10 GB)
# "mouse" downloads GRCm39-2024-A  (~10 GB)
SPECIES = "human"

# --- Samples ------------------------------------------------------------------
# Add one dictionary per sample.
#
# sample_id      : A short unique name. Used as the output folder name.
# sra_accession  : NCBI SRA run accession (e.g. "SRR12345678").
#                  Set to "" if your FASTQ files are already on disk.
# sample_prefix  : Must match the FASTQ filename prefix exactly.
#                  Cell Ranger expects:  <prefix>_S1_L001_R1_001.fastq.gz
#                  Usually the same as sample_id.

SAMPLES = [
    {
        "sample_id":      "pbmc_ctrl",
        "sra_accession":  "SRR12345678",
        "sample_prefix":  "pbmc_ctrl",
    },
    {
        "sample_id":      "pbmc_treat",
        "sra_accession":  "SRR12345679",
        "sample_prefix":  "pbmc_treat",
    },
]

# --- Cell Ranger parameters ---------------------------------------------------

# Number of CPU cores to use. More cores = faster. Check your machine first.
N_CORES = 16

# RAM in GB. Cell Ranger needs a lot — 64 GB minimum, 128 GB recommended.
MEM_GB = 128

# Expected number of cells per sample.
# Set to 0 to let Cell Ranger decide automatically (usually fine).
EXPECT_CELLS = 5000

# 10x chemistry version.
# "auto" works for most datasets. Options: SC3Pv2, SC3Pv3, SC3Pv3.1, SC5P-PE
CHEMISTRY = "auto"

# --- Optional: aggregate multiple samples -------------------------------------

# Set to True to run `cellranger aggr` after all samples finish.
# This merges samples but is NOT a substitute for Harmony batch correction
# (which is done in notebook 02). Usually you can leave this False.
RUN_AGGR = False
AGGR_ID = "combined"          # Name for the aggregated output folder
AGGR_NORMALIZE = "mapped"     # "mapped" normalises by sequencing depth; "none" skips

# --- Number of threads for download and QC steps ------------------------------
SRA_THREADS = 8      # fasterq-dump download threads
FASTQC_THREADS = 8   # FastQC threads


# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

def print_header(title: str) -> None:
    """Print a clearly visible section header."""
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def print_step(message: str) -> None:
    """Print a progress message with a bullet."""
    print(f"  >> {message}")


def print_ok(message: str) -> None:
    """Print a success message."""
    print(f"  [OK]  {message}")


def print_skip(message: str) -> None:
    """Print a skip message (step was already done)."""
    print(f"  [--]  {message}")


def print_error(message: str) -> None:
    """Print an error message."""
    print(f"  [!!]  {message}", file=sys.stderr)


def run_command(cmd: list, cwd: Path = None) -> None:
    """
    Run a shell command and stream its output to the terminal in real time.

    Parameters
    ----------
    cmd : list
        Command and arguments, e.g. ["cellranger", "count", "--id=s1", ...]
    cwd : Path, optional
        Directory to run the command from. Defaults to the current directory.

    Raises
    ------
    SystemExit
        If the command fails (non-zero exit code).
    """
    print_step("Running: " + " ".join(str(c) for c in cmd))
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        print_error(f"Command failed with exit code {result.returncode}.")
        print_error(f"Failed command: {' '.join(str(c) for c in cmd)}")
        sys.exit(1)


# ==============================================================================
#  STEP 1 — ENVIRONMENT CHECK
# ==============================================================================

def check_environment() -> None:
    """
    Verify that all required command-line tools are available.

    Cell Ranger, FastQC, and MultiQC are required.
    SRA Toolkit tools (fasterq-dump, prefetch) are required only if you are
    downloading data from NCBI SRA.
    """
    print_header("Step 1 — Environment Check")

    # tool name → is it required to run the pipeline?
    tools = {
        "cellranger":   True,
        "fastqc":       True,
        "multiqc":      True,
        "fasterq-dump": True,   # needed for SRA download
        "gzip":         True,
        "wget":         True,
        "samtools":     False,  # optional but useful
    }

    all_ok = True
    for tool, required in tools.items():
        found = shutil.which(tool) is not None
        status = "[OK]" if found else ("[!!] MISSING" if required else "[--] not found (optional)")
        print(f"  {status:20s}  {tool}")
        if required and not found:
            all_ok = False

    if not all_ok:
        print()
        print_error("One or more required tools are missing.")
        print_error("Install them and re-run this script.")
        print_error("  Cell Ranger  : https://www.10xgenomics.com/support/software/cell-ranger/downloads")
        print_error("  Others       : conda install -c bioconda sra-tools fastqc samtools")
        print_error("  MultiQC      : pip install multiqc")
        sys.exit(1)

    print_ok("All required tools found.")


# ==============================================================================
#  STEP 2 — REFERENCE GENOME
# ==============================================================================

# URLs for 10x Genomics pre-built reference packages
REFERENCE_URLS = {
    "human": (
        "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz",
        "refdata-gex-GRCh38-2024-A.tar.gz",
    ),
    "mouse": (
        "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz",
        "refdata-gex-GRCm39-2024-A.tar.gz",
    ),
}


def download_reference() -> None:
    """
    Download the Cell Ranger reference genome if it is not already on disk.

    The reference is ~10 GB and only needs to be downloaded once.
    It is used by Cell Ranger to align reads during `cellranger count`.
    """
    print_header("Step 2 — Reference Genome")

    if REFERENCE_PATH.exists():
        print_skip(f"Reference already exists at: {REFERENCE_PATH}")
        print_skip("Skipping download.")
        return

    if SPECIES not in REFERENCE_URLS:
        print_error(f"Unknown species '{SPECIES}'. Choose 'human' or 'mouse'.")
        sys.exit(1)

    url, tar_filename = REFERENCE_URLS[SPECIES]
    tar_path = REFERENCE_PATH.parent / tar_filename

    print_step(f"Downloading {SPECIES} reference genome (~10 GB) ...")
    print_step(f"Destination: {REFERENCE_PATH.parent}")

    REFERENCE_PATH.parent.mkdir(parents=True, exist_ok=True)

    # Download the compressed archive
    run_command(["wget", "-O", str(tar_path), url])

    # Extract it
    print_step("Extracting archive ...")
    run_command(["tar", "-xzvf", str(tar_path), "-C", str(REFERENCE_PATH.parent)])

    # Remove the .tar.gz to save disk space
    tar_path.unlink()
    print_ok(f"Reference genome ready at: {REFERENCE_PATH}")


# ==============================================================================
#  STEP 3 — FASTQ DOWNLOAD FROM SRA
# ==============================================================================

def download_fastqs() -> None:
    """
    Download FASTQ files from NCBI SRA for each sample using fasterq-dump.

    What happens for each sample:
      1. fasterq-dump downloads and splits the reads into _1.fastq and _2.fastq
      2. gzip compresses them to save disk space
      3. Files are renamed to match the 10x Cell Ranger naming convention:
            <prefix>_S1_L001_R1_001.fastq.gz   (barcode + UMI read)
            <prefix>_S1_L001_R2_001.fastq.gz   (cDNA read)

    If a sample has no SRA accession, or if the files already exist, it is skipped.
    """
    print_header("Step 3 — FASTQ Download from SRA")

    FASTQ_DIR.mkdir(parents=True, exist_ok=True)

    for sample in SAMPLES:
        sample_id  = sample["sample_id"]
        accession  = sample.get("sra_accession", "").strip()
        prefix     = sample.get("sample_prefix", sample_id).strip()

        # If there is no SRA accession, the user has their own FASTQ files
        if not accession:
            print_skip(f"{sample_id}: no SRA accession — assuming FASTQ files are already on disk")
            continue

        # Check if already downloaded (idempotent — safe to re-run)
        r1 = FASTQ_DIR / f"{prefix}_S1_L001_R1_001.fastq.gz"
        r2 = FASTQ_DIR / f"{prefix}_S1_L001_R2_001.fastq.gz"
        if r1.exists() and r2.exists():
            print_skip(f"{sample_id}: FASTQ files already present — skipping download")
            continue

        print_step(f"{sample_id}: downloading {accession} from NCBI SRA ...")

        # Download and split into separate read files
        run_command([
            "fasterq-dump",
            "--split-files",
            "--threads", str(SRA_THREADS),
            "--outdir", str(FASTQ_DIR),
            accession,
        ])

        # Compress and rename to Cell Ranger naming convention
        # SRA produces:  SRR12345678_1.fastq  and  SRR12345678_2.fastq
        # Cell Ranger needs:  prefix_S1_L001_R1_001.fastq.gz  etc.
        for read_number, cellranger_tag in [("1", "R1"), ("2", "R2")]:
            raw_file = FASTQ_DIR / f"{accession}_{read_number}.fastq"
            final_name = FASTQ_DIR / f"{prefix}_S1_L001_{cellranger_tag}_001.fastq.gz"

            if not raw_file.exists():
                print_error(f"Expected file not found after download: {raw_file}")
                sys.exit(1)

            # gzip compresses in-place and adds .gz extension
            run_command(["gzip", str(raw_file)])
            (FASTQ_DIR / f"{accession}_{read_number}.fastq.gz").rename(final_name)

        print_ok(f"{sample_id}: saved to {FASTQ_DIR}/{prefix}_S1_L001_R[12]_001.fastq.gz")


# ==============================================================================
#  STEP 4 — RAW READ QUALITY CONTROL
# ==============================================================================

def run_fastqc() -> None:
    """
    Run FastQC on every FASTQ file, then combine results with MultiQC.

    FastQC checks:
      - Per-base quality scores  (target: Q30 > 90% for R1 barcodes, > 85% for R2 cDNA)
      - GC content               (should match the expected GC% for your organism)
      - Adapter contamination    (< 5% is acceptable)
      - Sequence duplication     (scRNA-seq data naturally has high duplication from UMIs)

    MultiQC aggregates all individual FastQC reports into one HTML file.
    Open data/raw/fastqc_reports/multiqc_report.html in your browser to view results.
    """
    print_header("Step 4 — Raw Read QC (FastQC + MultiQC)")

    fastq_files = sorted(FASTQ_DIR.glob("*.fastq.gz"))

    if not fastq_files:
        print_error(f"No .fastq.gz files found in {FASTQ_DIR}")
        print_error("Run Step 3 first, or check that your FASTQ files are in the correct folder.")
        sys.exit(1)

    print_step(f"Found {len(fastq_files)} FASTQ file(s) to analyse.")

    qc_dir = FASTQ_DIR.parent / "fastqc_reports"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Run FastQC on every file
    run_command(
        ["fastqc", "--outdir", str(qc_dir), "--threads", str(FASTQC_THREADS)]
        + [str(f) for f in fastq_files]
    )

    # Combine all FastQC reports into one MultiQC report
    run_command([
        "multiqc", str(qc_dir),
        "--outdir", str(qc_dir),
        "--filename", "multiqc_report",
        "--force",   # overwrite if the report already exists
    ])

    report = qc_dir / "multiqc_report.html"
    print_ok(f"QC complete. Open this file in your browser to review:")
    print_ok(f"  {report}")


# ==============================================================================
#  STEP 5 — CELL RANGER COUNT
# ==============================================================================

def run_cellranger_count() -> None:
    """
    Run `cellranger count` for each sample to align reads and build count matrices.

    What Cell Ranger does internally:
      1. Aligns R2 (cDNA) reads to the reference genome using STAR
      2. Corrects cell barcodes against the 10x barcode whitelist
      3. Deduplicates UMIs (removes PCR duplicates)
      4. Counts how many unique RNA molecules (UMIs) came from each gene in each cell
      5. Writes results as a sparse matrix (genes × cells)

    Output per sample (inside data/raw/<sample_id>/outs/):
      filtered_feature_bc_matrix/  ← use this in notebook 01
      raw_feature_bc_matrix/
      web_summary.html             ← QC report; open in browser
      metrics_summary.csv          ← key QC numbers

    Runtime: roughly 2–6 hours per sample on 16 cores with 128 GB RAM.
    """
    print_header("Step 5 — Cell Ranger Count")

    CELLRANGER_OUT_DIR.mkdir(parents=True, exist_ok=True)

    for sample in SAMPLES:
        sample_id = sample["sample_id"]
        prefix    = sample.get("sample_prefix", sample_id).strip()

        # Check if this sample was already processed (idempotent)
        output_matrix = CELLRANGER_OUT_DIR / sample_id / "outs" / "filtered_feature_bc_matrix"
        if output_matrix.exists():
            print_skip(f"{sample_id}: output already exists — skipping")
            print_skip(f"  Delete {CELLRANGER_OUT_DIR / sample_id} to reprocess this sample")
            continue

        print_step(f"{sample_id}: starting Cell Ranger count ...")

        cmd = [
            "cellranger", "count",
            f"--id={sample_id}",
            f"--transcriptome={REFERENCE_PATH}",
            f"--fastqs={FASTQ_DIR}",
            f"--sample={prefix}",             # must match the FASTQ filename prefix
            f"--localcores={N_CORES}",
            f"--localmem={MEM_GB}",
        ]

        # Add optional parameters only if the user set them
        if CHEMISTRY != "auto":
            cmd.append(f"--chemistry={CHEMISTRY}")
        if EXPECT_CELLS > 0:
            cmd.append(f"--expect-cells={EXPECT_CELLS}")

        # Cell Ranger creates the output folder inside whatever directory it is
        # run from, so we run it from CELLRANGER_OUT_DIR
        run_command(cmd, cwd=CELLRANGER_OUT_DIR)
        print_ok(f"{sample_id}: count matrix written to {output_matrix}")


# ==============================================================================
#  STEP 6 — OPTIONAL MULTI-SAMPLE AGGREGATION
# ==============================================================================

def run_cellranger_aggr() -> None:
    """
    Merge all samples into a single count matrix using `cellranger aggr`.

    This step is OPTIONAL. It normalises sequencing depth across samples so that
    no single sample dominates the merged matrix. However, it does NOT correct
    for batch effects — use Harmony or scVI in notebook 02 for that.

    Tip: Most analyses skip this and load samples individually in notebook 01.
    """
    print_header("Step 6 — cellranger aggr (multi-sample merge)")

    if not RUN_AGGR:
        print_skip("RUN_AGGR = False — skipping. Set it to True in the configuration to enable.")
        return

    # Build the CSV file that cellranger aggr needs.
    # It lists each sample and the path to its molecule_info.h5 file.
    aggr_rows = []
    for sample in SAMPLES:
        sample_id = sample["sample_id"]
        h5_file = CELLRANGER_OUT_DIR / sample_id / "outs" / "molecule_info.h5"
        if h5_file.exists():
            aggr_rows.append({"sample_id": sample_id, "molecule_h5": str(h5_file)})
        else:
            print_error(f"{sample_id}: molecule_info.h5 not found at {h5_file}")
            print_error("Run Step 5 (Cell Ranger count) for all samples first.")
            sys.exit(1)

    csv_path = CELLRANGER_OUT_DIR / "aggr_samples.csv"
    pd.DataFrame(aggr_rows).to_csv(csv_path, index=False)
    print_step(f"Aggregation CSV written: {csv_path}")

    run_command([
        "cellranger", "aggr",
        f"--id={AGGR_ID}",
        f"--csv={csv_path}",
        f"--normalize={AGGR_NORMALIZE}",
    ], cwd=CELLRANGER_OUT_DIR)

    out = CELLRANGER_OUT_DIR / AGGR_ID / "outs"
    print_ok(f"Aggregation complete: {out}")


# ==============================================================================
#  STEP 7 — QC METRICS SUMMARY
# ==============================================================================

def print_metrics_summary() -> None:
    """
    Read Cell Ranger's metrics_summary.csv for each sample and print a table.

    Key metrics to check before moving to notebook 01:

    Metric                      What to look for
    ─────────────────────────   ────────────────────────────────────────────
    Estimated number of cells   Within ±20% of expected
    Mean reads per cell         20,000–50,000 (more = better coverage)
    Median genes per cell       1,000–5,000   (< 500 suggests problems)
    Sequencing saturation       70–85%        (> 85% = diminishing returns)
    Q30 bases in barcode        > 90%
    Q30 bases in RNA read       > 85%
    """
    print_header("Step 7 — QC Metrics Summary")

    all_metrics = []
    for sample in SAMPLES:
        sample_id = sample["sample_id"]
        csv_path = CELLRANGER_OUT_DIR / sample_id / "outs" / "metrics_summary.csv"
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            df.insert(0, "Sample", sample_id)
            all_metrics.append(df)
        else:
            print_skip(f"{sample_id}: metrics_summary.csv not found — run Step 5 first")

    if not all_metrics:
        print_error("No metrics files found. Run Cell Ranger count first (Step 5).")
        return

    combined = pd.concat(all_metrics, ignore_index=True)

    # Print the table in a readable format
    print()
    with pd.option_context("display.max_columns", None, "display.width", 120):
        print(combined.to_string(index=False))


# ==============================================================================
#  OUTPUT SUMMARY
# ==============================================================================

def print_output_summary() -> None:
    """
    Print the output paths for each sample so you know what to pass to notebook 01.
    """
    print_header("Output Summary — What to pass to notebook 01")

    print()
    print(f"  {'Sample':<20}  {'filtered_matrix_path'}")
    print(f"  {'-'*20}  {'-'*50}")

    for sample in SAMPLES:
        sample_id = sample["sample_id"]
        matrix_path = CELLRANGER_OUT_DIR / sample_id / "outs" / "filtered_feature_bc_matrix"
        if matrix_path.exists():
            print(f"  {sample_id:<20}  {matrix_path}")
        else:
            print(f"  {sample_id:<20}  (not found — Cell Ranger not run yet)")

    print()
    print("  In notebook 01, set 'Input path' to the path shown above for your sample.")
    print("  Run notebook 01 once per sample, then merge them in notebook 02.")
    print()
    print("  Also check the web_summary.html files in each sample's outs/ folder.")
    print("  Open them in your browser for a detailed interactive QC report.")
    print()
    print("  Files you should keep after this step:")
    print("    filtered_feature_bc_matrix/   <- input to notebook 01 (required)")
    print("    raw_feature_bc_matrix/        <- useful for ambient RNA analysis")
    print("    molecule_info.h5              <- needed if you run cellranger aggr")
    print("    web_summary.html              <- keep for publication records")
    print("    metrics_summary.csv           <- keep for publication records")
    print("    possorted_genome_bam.bam      <- optional (~20-50 GB); needed for RNA velocity")
    print("    SC_RNA_COUNTER_CS/            <- safe to delete (intermediate files)")


# ==============================================================================
#  MAIN — Runs all steps in order
# ==============================================================================

def main() -> None:
    print()
    print("=" * 70)
    print("  scRNA-seq Preprocessing: FASTQ → Count Matrix")
    print("=" * 70)
    print(f"  Species  : {SPECIES}")
    print(f"  Samples  : {len(SAMPLES)}")
    print(f"  FASTQ dir: {FASTQ_DIR}")
    print(f"  Output   : {CELLRANGER_OUT_DIR}")
    print(f"  Reference: {REFERENCE_PATH}")
    print(f"  Cores    : {N_CORES}   RAM: {MEM_GB} GB")

    check_environment()      # Step 1: make sure tools are installed
    download_reference()     # Step 2: get the reference genome (skips if already done)
    download_fastqs()        # Step 3: download from SRA (skips if already done)
    run_fastqc()             # Step 4: check raw read quality
    run_cellranger_count()   # Step 5: align and count (skips completed samples)
    run_cellranger_aggr()    # Step 6: merge samples (skipped if RUN_AGGR = False)
    print_metrics_summary()  # Step 7: show QC numbers
    print_output_summary()   # Show paths to pass to notebook 01

    print()
    print("  Done. Open notebook 01 (01_qc_filtering.py) for the next step.")
    print()


if __name__ == "__main__":
    main()
