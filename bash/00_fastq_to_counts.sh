#!/usr/bin/env bash
# 00_fastq_to_counts.sh
# ======================
# FASTQ -> Count Matrix preprocessing pipeline for 10x Chromium scRNA-seq data.
#
# WHAT THIS SCRIPT DOES (runs top-to-bottom in this order):
#   Step 1 — Check that all required tools are installed
#   Step 2 — Download the Cell Ranger reference genome (only needed once, ~10 GB)
#   Step 3 — Download FASTQ files from NCBI SRA (skip if files already exist)
#   Step 4 — Run FastQC + MultiQC to assess raw read quality
#   Step 5 — Run Cell Ranger count for each sample
#   Step 6 — (Optional) Merge all samples with cellranger aggr
#   Step 7 — Print QC metrics summary
#
# HOW TO USE:
#   1. Edit the CONFIGURATION section below.
#   2. Make executable (once):  chmod +x 00_fastq_to_counts.sh
#   3. Run:                     bash 00_fastq_to_counts.sh
#
# REQUIRED TOOLS (install once before running):
#   Cell Ranger >= 8.0  https://www.10xgenomics.com/support/software/cell-ranger/downloads
#   SRA Toolkit         pixi install -e genomics  (or: conda install -c bioconda sra-tools)
#   FastQC              pixi install -e genomics  (or: conda install -c bioconda fastqc)
#   MultiQC             pixi install -e genomics  (or: conda install -c bioconda multiqc)
#   samtools            pixi install -e genomics  (or: conda install -c bioconda samtools)
#
# OUTPUT:
#   data/raw/<sample_id>/outs/filtered_feature_bc_matrix/
#   data/raw/fastqc_reports/multiqc_report.html

set -euo pipefail


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

# --- Paths --------------------------------------------------------------------

# Where the Cell Ranger pre-built reference is (or will be downloaded to).
# Only downloaded once and reusable across projects.
REFERENCE_PATH="${HOME}/references/cellranger/refdata-gex-GRCh38-2024-A"

# Folder where FASTQ files will be stored / are already present.
FASTQ_DIR="data/raw/fastqs"

# Root folder for Cell Ranger output. Each sample gets <OUTPUT_DIR>/<sample_id>/outs/
CELLRANGER_OUT_DIR="data/raw"

# --- Reference genome ---------------------------------------------------------

# "human" downloads refdata-gex-GRCh38-2024-A  (~10 GB)
# "mouse" downloads refdata-gex-GRCm39-2024-A  (~10 GB)
SPECIES="human"

# --- Samples ------------------------------------------------------------------
# Three parallel arrays — one entry per sample (same index = same sample).
#
#   SAMPLE_IDS       short unique name; becomes the output folder name
#   SRA_ACCESSIONS   NCBI SRA run accession (e.g. SRR12345678)
#                    Set to "" if FASTQs are already on disk
#   SAMPLE_PREFIXES  must match the FASTQ filename prefix exactly
#                    Cell Ranger expects: <prefix>_S1_L001_R[12]_001.fastq.gz

SAMPLE_IDS=(        "pbmc_ctrl"    "pbmc_treat"   )
SRA_ACCESSIONS=(    "SRR12345678"  "SRR12345679"  )
SAMPLE_PREFIXES=(   "pbmc_ctrl"    "pbmc_treat"   )

# --- Cell Ranger parameters ---------------------------------------------------

N_CORES=16           # CPU cores
MEM_GB=128           # RAM in GB (minimum 64, recommended 128)
EXPECT_CELLS=5000    # set to 0 to let Cell Ranger decide automatically
CHEMISTRY="auto"     # auto works for most datasets; options: SC3Pv2 SC3Pv3 SC3Pv3.1 SC5P-PE

# --- SRA / QC threads ---------------------------------------------------------

SRA_THREADS=8
FASTQC_THREADS=8

# --- Optional aggregation -----------------------------------------------------

RUN_AGGR=false       # set to true to run cellranger aggr after all samples finish
AGGR_ID="combined"
AGGR_NORMALIZE="mapped"   # "mapped" normalises by sequencing depth; "none" skips


# ==============================================================================
#  REFERENCE GENOME URLS
# ==============================================================================

HUMAN_REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
MOUSE_REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"


# ==============================================================================
#  HELPERS
# ==============================================================================

header() { echo; echo "======================================================================"; echo "  $1"; echo "======================================================================"; }
ok()     { echo "  [OK]  $1"; }
skip()   { echo "  [--]  $1"; }
info()   { echo "  >>    $1"; }
err()    { echo "  [!!]  $1" >&2; }


# ==============================================================================
#  STEP 1 — ENVIRONMENT CHECK
# ==============================================================================

check_environment() {
    header "Step 1 — Environment Check"

    local all_ok=true
    local tools_required=(cellranger fastqc multiqc fasterq-dump gzip wget)
    local tools_optional=(samtools)

    for tool in "${tools_required[@]}"; do
        if command -v "$tool" &>/dev/null; then
            ok "$tool"
        else
            err "MISSING (required): $tool"
            all_ok=false
        fi
    done

    for tool in "${tools_optional[@]}"; do
        if command -v "$tool" &>/dev/null; then
            ok "$tool"
        else
            skip "$tool  (optional — not required)"
        fi
    done

    if [[ "$all_ok" == false ]]; then
        err "One or more required tools are missing."
        err "  Cell Ranger : https://www.10xgenomics.com/support/software/cell-ranger/downloads"
        err "  Others      : pixi install -e genomics"
        exit 1
    fi
    ok "All required tools found."
}


# ==============================================================================
#  STEP 2 — REFERENCE GENOME
# ==============================================================================

download_reference() {
    header "Step 2 — Reference Genome"

    if [[ -d "$REFERENCE_PATH" ]]; then
        skip "Reference already exists at: ${REFERENCE_PATH}"
        skip "Skipping download."
        return
    fi

    local url tar_name ref_parent
    ref_parent="$(dirname "$REFERENCE_PATH")"

    case "$SPECIES" in
        human) url="$HUMAN_REF_URL" ; tar_name="refdata-gex-GRCh38-2024-A.tar.gz" ;;
        mouse) url="$MOUSE_REF_URL" ; tar_name="refdata-gex-GRCm39-2024-A.tar.gz" ;;
        *)
            err "Unknown SPECIES '${SPECIES}'. Choose 'human' or 'mouse'."
            exit 1 ;;
    esac

    local tar_path="${ref_parent}/${tar_name}"
    mkdir -p "$ref_parent"

    info "Downloading ${SPECIES} reference genome (~10 GB) ..."
    info "Destination: ${ref_parent}"
    wget -O "$tar_path" "$url"

    info "Extracting archive ..."
    tar -xzf "$tar_path" -C "$ref_parent"
    rm -f "$tar_path"

    ok "Reference genome ready at: ${REFERENCE_PATH}"
}


# ==============================================================================
#  STEP 3 — FASTQ DOWNLOAD FROM SRA
# ==============================================================================

download_fastqs() {
    header "Step 3 — FASTQ Download from SRA"

    mkdir -p "$FASTQ_DIR"
    local n="${#SAMPLE_IDS[@]}"

    for (( i=0; i<n; i++ )); do
        local sample_id="${SAMPLE_IDS[$i]}"
        local accession="${SRA_ACCESSIONS[$i]}"
        local prefix="${SAMPLE_PREFIXES[$i]}"

        if [[ -z "$accession" ]]; then
            skip "${sample_id}: no SRA accession — assuming FASTQs are already on disk"
            continue
        fi

        local r1="${FASTQ_DIR}/${prefix}_S1_L001_R1_001.fastq.gz"
        local r2="${FASTQ_DIR}/${prefix}_S1_L001_R2_001.fastq.gz"

        if [[ -f "$r1" && -f "$r2" ]]; then
            skip "${sample_id}: FASTQ files already present — skipping download"
            continue
        fi

        info "${sample_id}: downloading ${accession} from NCBI SRA ..."
        fasterq-dump \
            --split-files \
            --threads "$SRA_THREADS" \
            --outdir  "$FASTQ_DIR" \
            "$accession"

        # Compress and rename to Cell Ranger naming convention
        # fasterq-dump produces:  SRR12345678_1.fastq  SRR12345678_2.fastq
        # Cell Ranger needs:      prefix_S1_L001_R1_001.fastq.gz  etc.
        for read_num in 1 2; do
            local raw="${FASTQ_DIR}/${accession}_${read_num}.fastq"
            local read_tag="R${read_num}"
            local final="${FASTQ_DIR}/${prefix}_S1_L001_${read_tag}_001.fastq.gz"

            if [[ ! -f "$raw" ]]; then
                err "Expected file not found after download: ${raw}"
                exit 1
            fi

            gzip "$raw"
            mv "${raw}.gz" "$final"
        done

        ok "${sample_id}: ${FASTQ_DIR}/${prefix}_S1_L001_R[12]_001.fastq.gz"
    done
}


# ==============================================================================
#  STEP 4 — RAW READ QC (FastQC + MultiQC)
# ==============================================================================

run_fastqc() {
    header "Step 4 — Raw Read QC (FastQC + MultiQC)"

    mapfile -t fastq_files < <(find "$FASTQ_DIR" -name "*.fastq.gz" | sort)

    if [[ "${#fastq_files[@]}" -eq 0 ]]; then
        err "No .fastq.gz files found in ${FASTQ_DIR}"
        err "Run Step 3 first, or verify your FASTQ files are in the correct folder."
        exit 1
    fi

    info "Found ${#fastq_files[@]} FASTQ file(s)."

    local qc_dir="${CELLRANGER_OUT_DIR}/fastqc_reports"
    mkdir -p "$qc_dir"

    fastqc \
        --outdir  "$qc_dir" \
        --threads "$FASTQC_THREADS" \
        "${fastq_files[@]}"

    multiqc "$qc_dir" \
        --outdir   "$qc_dir" \
        --filename multiqc_report \
        --force

    ok "QC complete. Open in browser: ${qc_dir}/multiqc_report.html"
}


# ==============================================================================
#  STEP 5 — CELL RANGER COUNT
# ==============================================================================

run_cellranger_count() {
    header "Step 5 — Cell Ranger Count"

    if [[ ! -d "$REFERENCE_PATH" ]]; then
        err "Reference not found: ${REFERENCE_PATH}"
        err "Run Step 2 first."
        exit 1
    fi

    mkdir -p "$CELLRANGER_OUT_DIR"
    local n="${#SAMPLE_IDS[@]}"

    for (( i=0; i<n; i++ )); do
        local sample_id="${SAMPLE_IDS[$i]}"
        local prefix="${SAMPLE_PREFIXES[$i]}"
        local output_matrix="${CELLRANGER_OUT_DIR}/${sample_id}/outs/filtered_feature_bc_matrix"

        if [[ -d "$output_matrix" ]]; then
            skip "${sample_id}: output already exists — skipping"
            skip "  Delete ${CELLRANGER_OUT_DIR}/${sample_id} to reprocess"
            continue
        fi

        info "${sample_id}: starting Cell Ranger count ..."

        local cmd=(
            cellranger count
            "--id=${sample_id}"
            "--transcriptome=${REFERENCE_PATH}"
            "--fastqs=$(realpath "$FASTQ_DIR")"
            "--sample=${prefix}"
            "--localcores=${N_CORES}"
            "--localmem=${MEM_GB}"
        )

        [[ "$CHEMISTRY"    != "auto" ]] && cmd+=( "--chemistry=${CHEMISTRY}" )
        [[ "$EXPECT_CELLS" -gt 0    ]] && cmd+=( "--expect-cells=${EXPECT_CELLS}" )

        # Cell Ranger writes output relative to CWD, so run from the output dir
        (cd "$CELLRANGER_OUT_DIR" && "${cmd[@]}")

        ok "${sample_id}: count matrix written to ${output_matrix}"
    done
}


# ==============================================================================
#  STEP 6 — CELL RANGER AGGR (OPTIONAL)
# ==============================================================================

run_cellranger_aggr() {
    header "Step 6 — cellranger aggr (multi-sample merge)"

    if [[ "$RUN_AGGR" != true ]]; then
        skip "RUN_AGGR=false — skipping. Set RUN_AGGR=true to enable."
        return
    fi

    local csv_path="${CELLRANGER_OUT_DIR}/aggr_samples.csv"
    local n="${#SAMPLE_IDS[@]}"

    # Build the aggr CSV: header + one row per sample
    echo "sample_id,molecule_h5" > "$csv_path"

    for (( i=0; i<n; i++ )); do
        local sample_id="${SAMPLE_IDS[$i]}"
        local h5="${CELLRANGER_OUT_DIR}/${sample_id}/outs/molecule_info.h5"

        if [[ ! -f "$h5" ]]; then
            err "${sample_id}: molecule_info.h5 not found at ${h5}"
            err "Run Step 5 for all samples first."
            exit 1
        fi
        echo "${sample_id},${h5}" >> "$csv_path"
    done

    info "Aggregation CSV written: ${csv_path}"

    (cd "$CELLRANGER_OUT_DIR" && cellranger aggr \
        "--id=${AGGR_ID}"       \
        "--csv=${csv_path}"     \
        "--normalize=${AGGR_NORMALIZE}")

    ok "Aggregation complete: ${CELLRANGER_OUT_DIR}/${AGGR_ID}/outs/"
}


# ==============================================================================
#  STEP 7 — QC METRICS SUMMARY
# ==============================================================================
#
#  Key metrics to check before moving to script 01:
#
#  Metric                       What to look for
#  ---------------------------  -------------------------------------------------
#  Estimated number of cells    Within ±20% of expected
#  Mean reads per cell          20,000–50,000  (more = better coverage)
#  Median genes per cell        1,000–5,000    (< 500 suggests problems)
#  Sequencing saturation        70–85%         (> 85% = diminishing returns)
#  Q30 bases in barcode         > 90%
#  Q30 bases in RNA read        > 85%

print_metrics_summary() {
    header "Step 7 — QC Metrics Summary"

    local n="${#SAMPLE_IDS[@]}"
    local found=false

    for (( i=0; i<n; i++ )); do
        local sample_id="${SAMPLE_IDS[$i]}"
        local csv="${CELLRANGER_OUT_DIR}/${sample_id}/outs/metrics_summary.csv"

        if [[ -f "$csv" ]]; then
            echo
            echo "  Sample: ${sample_id}"
            echo "  ----------------------------------------"
            # Print each metric on its own line (CSV has header + one data row)
            paste -d '=' \
                <(head -1 "$csv" | tr ',' '\n') \
                <(tail -1 "$csv" | tr ',' '\n') \
            | sed 's/^/    /'
            found=true
        else
            skip "${sample_id}: metrics_summary.csv not found — run Step 5 first"
        fi
    done

    if [[ "$found" == false ]]; then
        err "No metrics files found. Run Cell Ranger count first (Step 5)."
    fi
}


# ==============================================================================
#  OUTPUT SUMMARY
# ==============================================================================

print_output_summary() {
    header "Output Summary — Paths to pass to script 01"

    echo
    printf "  %-22s  %s\n" "Sample" "filtered_matrix_path"
    printf "  %-22s  %s\n" "----------------------" "--------------------------------------------------"

    local n="${#SAMPLE_IDS[@]}"
    for (( i=0; i<n; i++ )); do
        local sample_id="${SAMPLE_IDS[$i]}"
        local matrix="${CELLRANGER_OUT_DIR}/${sample_id}/outs/filtered_feature_bc_matrix"
        if [[ -d "$matrix" ]]; then
            printf "  %-22s  %s\n" "$sample_id" "$matrix"
        else
            printf "  %-22s  (not found — Cell Ranger not run yet)\n" "$sample_id"
        fi
    done

    echo
    echo "  Run script 01 once per sample, then merge them in script 02."
    echo
    echo "  Files to keep after this step:"
    echo "    filtered_feature_bc_matrix/   <- input to script 01 (required)"
    echo "    raw_feature_bc_matrix/        <- useful for ambient RNA analysis"
    echo "    molecule_info.h5              <- needed if you run cellranger aggr"
    echo "    web_summary.html              <- keep for publication records"
    echo "    metrics_summary.csv           <- keep for publication records"
    echo "    possorted_genome_bam.bam      <- optional; needed for RNA velocity"
    echo "    SC_RNA_COUNTER_CS/            <- safe to delete (intermediate files)"
    echo
}


# ==============================================================================
#  MAIN — Runs all steps in order
# ==============================================================================

main() {
    echo
    echo "======================================================================"
    echo "  scRNA-seq Preprocessing: FASTQ -> Count Matrix"
    echo "======================================================================"
    echo "  Species    : ${SPECIES}"
    echo "  Samples    : ${#SAMPLE_IDS[@]}"
    echo "  FASTQ dir  : ${FASTQ_DIR}"
    echo "  Output dir : ${CELLRANGER_OUT_DIR}"
    echo "  Reference  : ${REFERENCE_PATH}"
    echo "  Cores      : ${N_CORES}   RAM: ${MEM_GB} GB"

    check_environment        # Step 1: verify tools
    download_reference       # Step 2: get reference genome (skips if already done)
    download_fastqs          # Step 3: download from SRA  (skips if already done)
    run_fastqc               # Step 4: assess raw read quality
    run_cellranger_count     # Step 5: align and count    (skips completed samples)
    run_cellranger_aggr      # Step 6: merge samples      (skipped if RUN_AGGR=false)
    print_metrics_summary    # Step 7: show QC numbers
    print_output_summary     # Show paths to pass to script 01

    echo
    echo "  Done. Open script 01 (01_qc_filtering.R) for the next step."
    echo
}

main "$@"
