# 07_cnv_analysis.R
# ==================
# Copy number variation (CNV) inference using inferCNV.
#
# WHAT THIS SCRIPT DOES:
#   1. Load annotated tumour + normal data
#   2. Build the cell annotation file required by inferCNV
#   3. Run inferCNV to infer CNVs from expression patterns
#   4. Classify cells as aneuploid (tumour) or diploid (normal) by CNV score
#   5. Transfer CNV status back to the Seurat object and save
#
# WHEN TO USE:
#   Tumour datasets where you need to distinguish malignant cells from normal
#   stromal/immune cells, or to identify tumour subclones.
#
# IMPORTANT CAVEATS:
#   - inferCNV infers CNV from expression — ~80% concordance with DNA WGS.
#   - Always interpret in context of cell type annotations.
#   - A matched normal sample dramatically improves specificity.
#
# HOW TO USE:
#   1. Edit CONFIGURATION below — especially NORMAL_CELL_TYPES and TUMOUR_CELL_TYPES.
#   2. Run:  Rscript 07_cnv_analysis.R
#
# OUTPUT:
#   results/07_cnv/seurat_cnv.rds
#   results/07_cnv/infercnv_output/   (all inferCNV files)
#   results/07_cnv/plots/

suppressMessages({
  library(Seurat)
  library(infercnv)
  library(ggplot2)
})


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH    <- "results/03_annotation/seurat_annotated.rds"
CELL_TYPE_COL <- "cell_type"

# --- Reference (normal) and tumour populations --------------------------------
# Normal cells are used as the baseline for CNV inference.
# These should be non-malignant cells (immune, stromal, etc.).
NORMAL_CELL_TYPES <- c("T cell", "B cell", "NK cell")

# Cell types you suspect are tumour / malignant.
TUMOUR_CELL_TYPES <- c("Epithelial")

# --- Gene order file ----------------------------------------------------------
# inferCNV requires a gene ordering file (hg38 or mm10).
# Download from:
#   Human (hg38): https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt
#   Mouse (mm10): https://data.broadinstitute.org/Trinity/CTAT/cnv/mm10_gencode_vM25.txt
GENE_ORDER_FILE <- "data/raw/hg38_gencode_v27.txt"

# Reference genome. Used for inferCNV annotation.
# Options: "hg38", "hg19", "mm10", "mm9"
GENOME_VERSION <- "hg38"

# --- inferCNV parameters ------------------------------------------------------
# Number of genes per smoothing window.
WINDOW_SIZE <- 101L

# Noise threshold: expression differences smaller than this are treated as noise.
NOISE_FILTER <- 0.1

# CNV score threshold for aneuploid classification.
# Cells with variance of CNV profile > (median_normal + threshold) = Aneuploid.
CNV_SCORE_THRESHOLD <- 0.02

# Number of CPU threads.
N_THREADS <- 4L

# --- Output -------------------------------------------------------------------
OUTPUT_DIR  <- "results/07_cnv"
OUTPUT_PATH <- file.path(OUTPUT_DIR, "seurat_cnv.rds")
INFERCNV_DIR <- file.path(OUTPUT_DIR, "infercnv_output")


# ==============================================================================
#  STEP 1 — LOAD & VALIDATE
# ==============================================================================

load_data <- function() {
  if (!file.exists(INPUT_PATH)) stop("File not found: ", INPUT_PATH)
  seu <- readRDS(INPUT_PATH)

  all_types <- unique(as.character(seu@meta.data[[CELL_TYPE_COL]]))

  missing_normal <- setdiff(NORMAL_CELL_TYPES, all_types)
  missing_tumour <- setdiff(TUMOUR_CELL_TYPES, all_types)

  if (length(missing_normal) > 0L) {
    stop(sprintf("Normal cell types not found: %s\nAvailable: %s",
                 paste(missing_normal, collapse = ", "),
                 paste(all_types, collapse = ", ")))
  }
  if (length(missing_tumour) > 0L) {
    stop(sprintf("Tumour cell types not found: %s\nAvailable: %s",
                 paste(missing_tumour, collapse = ", "),
                 paste(all_types, collapse = ", ")))
  }
  if (!file.exists(GENE_ORDER_FILE)) {
    stop("Gene order file not found: ", GENE_ORDER_FILE,
         "\nDownload from https://data.broadinstitute.org/Trinity/CTAT/cnv/")
  }
  seu
}


# ==============================================================================
#  STEP 2 — BUILD CELL ANNOTATION FILE
# ==============================================================================

build_cell_annotations <- function(seu, out_dir) {
  all_types <- c(NORMAL_CELL_TYPES, TUMOUR_CELL_TYPES)
  mask      <- seu@meta.data[[CELL_TYPE_COL]] %in% all_types
  sub       <- seu[, mask]

  annot_df <- data.frame(
    cell_barcode = colnames(sub),
    cell_type    = as.character(sub@meta.data[[CELL_TYPE_COL]]),
    row.names    = NULL
  )

  annot_path <- file.path(out_dir, "cell_annotations.txt")
  write.table(annot_df, annot_path, sep = "\t",
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  list(sub = sub, annot_path = annot_path)
}


# ==============================================================================
#  STEP 3 — RUN inferCNV
# ==============================================================================

run_infercnv <- function(sub, annot_path) {
  dir.create(INFERCNV_DIR, recursive = TRUE, showWarnings = FALSE)

  # Use raw counts
  counts <- GetAssayData(sub, layer = "counts")

  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix    = counts,
    annotations_file     = annot_path,
    delim                = "\t",
    gene_order_file      = GENE_ORDER_FILE,
    ref_group_names      = NORMAL_CELL_TYPES
  )

  infercnv_obj <- run(
    infercnv_obj,
    cutoff               = NOISE_FILTER,
    out_dir              = INFERCNV_DIR,
    cluster_by_groups    = TRUE,
    plot_steps           = FALSE,
    denoise              = TRUE,
    HMM                  = FALSE,
    num_threads          = N_THREADS,
    window_length        = WINDOW_SIZE,
    output_format        = "pdf"
  )

  infercnv_obj
}


# ==============================================================================
#  STEP 4 — CNV SCORE & CLASSIFICATION
# ==============================================================================

classify_cnv <- function(seu, infercnv_obj) {
  # Extract the smoothed expression matrix (CNV signal)
  cnv_mat <- infercnv_obj@expr.data

  # CNV score = variance of CNV profile per cell (higher = more altered genome)
  shared_cells <- intersect(colnames(seu), colnames(cnv_mat))
  cnv_score    <- apply(cnv_mat[, shared_cells, drop = FALSE], 2, var)

  # Reference = normal cell types
  ref_mask     <- seu@meta.data[[CELL_TYPE_COL]] %in% NORMAL_CELL_TYPES
  ref_cells    <- intersect(colnames(seu)[ref_mask], names(cnv_score))
  median_ref   <- median(cnv_score[ref_cells], na.rm = TRUE)

  # Classify
  cnv_status                <- rep(NA_character_, ncol(seu))
  names(cnv_status)         <- colnames(seu)
  cnv_score_full            <- rep(NA_real_,      ncol(seu))
  names(cnv_score_full)     <- colnames(seu)

  cnv_score_full[names(cnv_score)] <- cnv_score
  cnv_status[names(cnv_score)]     <- ifelse(
    cnv_score > median_ref + CNV_SCORE_THRESHOLD,
    "Aneuploid", "Diploid"
  )

  seu$cnv_score  <- cnv_score_full
  seu$cnv_status <- cnv_status

  n_aneupl <- sum(cnv_status == "Aneuploid", na.rm = TRUE)
  n_total  <- sum(!is.na(cnv_status))
  message(sprintf("Aneuploid cells: %d / %d (%.1f%%)",
                  n_aneupl, n_total, 100 * n_aneupl / n_total))
  seu
}


# ==============================================================================
#  STEP 5 — VISUALISATION
# ==============================================================================

save_cnv_plots <- function(seu, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  if ("umap" %in% names(seu@reductions) && "cnv_status" %in% colnames(seu@meta.data)) {
    p_status <- DimPlot(
      seu, reduction = "umap", group.by = "cnv_status", pt.size = 0.4,
      cols = c("Aneuploid" = "#C44E52", "Diploid" = "#4C72B0")
    ) + ggtitle("CNV Status on UMAP") + theme_classic()
    ggsave(file.path(plot_dir, "01_umap_cnv_status.pdf"), p_status, width = 8, height = 7)

    p_score <- FeaturePlot(seu, features = "cnv_score", reduction = "umap",
                           pt.size = 0.4) +
      scale_colour_gradient(low = "grey90", high = "#C44E52") +
      ggtitle("CNV Score on UMAP") + theme_classic()
    ggsave(file.path(plot_dir, "02_umap_cnv_score.pdf"), p_score, width = 8, height = 7)
  }
}


# ==============================================================================
#  MAIN
# ==============================================================================

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(OUTPUT_DIR, "plots")

  seu                    <- load_data()
  result                 <- build_cell_annotations(seu, OUTPUT_DIR)
  infercnv_obj           <- run_infercnv(result$sub, result$annot_path)
  seu                    <- classify_cnv(seu, infercnv_obj)
  save_cnv_plots(seu, plot_dir)

  saveRDS(seu, OUTPUT_PATH)
  message("Saved -> ", OUTPUT_PATH)
}

main()
