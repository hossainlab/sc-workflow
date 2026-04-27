# 01_qc_filtering.R
# ==================
# Quality control and cell filtering for a single 10x scRNA-seq sample.
#
# WHAT THIS SCRIPT DOES:
#   1. Load Cell Ranger output (MTX directory or .h5 file)
#   2. Calculate per-cell QC metrics (UMI counts, gene counts, % mitochondrial)
#   3. Detect and remove doublets with scDblFinder
#   4. Filter cells and genes using the thresholds set below
#   5. Save a clean Seurat object to results/01_qc/
#
# HOW TO USE:
#   1. Edit the CONFIGURATION section below.
#   2. Run:  Rscript 01_qc_filtering.R
#   3. Run this script once per sample, then merge in script 02.
#
# OUTPUT:
#   results/01_qc/<sample_name>_qc.rds
#   results/01_qc/plots/

suppressMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

# Path to Cell Ranger output directory (filtered_feature_bc_matrix) or a .h5 file.
DATA_PATH <- "data/raw/sample1/outs/filtered_feature_bc_matrix"

# A short name for this sample. Used in output file names and plot titles.
SAMPLE_NAME <- "sample1"

# "human" or "mouse" — controls mitochondrial gene prefix (MT- vs mt-)
SPECIES <- "human"

# --- Cell filtering thresholds ------------------------------------------------
MIN_GENES    <- 500     # minimum genes detected per cell
MAX_GENES    <- 8000    # maximum genes (very high = likely doublet)
MIN_COUNTS   <- 1000    # minimum total UMI counts per cell
MAX_COUNTS   <- 50000   # maximum total UMI counts per cell
MAX_MT_PCT   <- 20      # maximum mitochondrial gene % (high = dying cell)

# Genes expressed in fewer than this many cells are removed.
MIN_CELLS_PER_GENE <- 3L

# --- Doublet detection --------------------------------------------------------
RUN_SCDBLFINDER    <- TRUE
EXPECTED_DBR       <- 0.06   # expected doublet rate (0.8% per 1000 cells loaded)

# --- Output -------------------------------------------------------------------
OUTPUT_DIR <- "results/01_qc"


# ==============================================================================
#  HELPERS
# ==============================================================================

mt_pattern <- function(species) {
  if (tolower(species) == "mouse") "^mt-" else "^MT-"
}

ribo_pattern <- function(species) {
  if (tolower(species) == "mouse") "^Rp[sl]" else "^RP[SL]"
}

hb_pattern <- function(species) {
  if (tolower(species) == "mouse") "^Hb[^p]" else "^HB[^P]"
}


# ==============================================================================
#  STEP 1 — LOAD DATA
# ==============================================================================

load_data <- function() {
  p <- DATA_PATH
  if (!file.exists(p) && !dir.exists(p)) {
    stop("Data path not found: ", p, "\nSet DATA_PATH in the configuration section.")
  }

  if (grepl("\\.h5$", p)) {
    counts <- Read10X_h5(p)
    if (is.list(counts)) counts <- counts[["Gene Expression"]]
  } else {
    counts <- Read10X(p)
  }

  seu <- CreateSeuratObject(
    counts  = counts,
    project = SAMPLE_NAME,
    min.cells = 0L,
    min.features = 0L
  )
  seu$sample <- SAMPLE_NAME
  message(sprintf("Loaded %d cells x %d genes", ncol(seu), nrow(seu)))
  seu
}


# ==============================================================================
#  STEP 2 — QC METRICS
# ==============================================================================

compute_qc_metrics <- function(seu) {
  seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = mt_pattern(SPECIES))
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = ribo_pattern(SPECIES))
  seu[["percent.hb"]]   <- PercentageFeatureSet(seu, pattern = hb_pattern(SPECIES))
  seu[["log10_counts"]] <- log10(seu$nCount_RNA + 1)
  seu[["log10_genes"]]  <- log10(seu$nFeature_RNA + 1)

  message(sprintf(
    "Median UMI/cell: %.0f  |  Median genes/cell: %.0f  |  Median MT%%: %.2f%%",
    median(seu$nCount_RNA),
    median(seu$nFeature_RNA),
    median(seu$percent.mt)
  ))
  seu
}

save_qc_plots <- function(seu, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  p1 <- VlnPlot(seu, features = "nCount_RNA",   pt.size = 0, ncol = 1) +
    ggtitle("UMI Counts") + NoLegend()
  p2 <- VlnPlot(seu, features = "nFeature_RNA", pt.size = 0, ncol = 1) +
    ggtitle("Genes Detected") + NoLegend()
  p3 <- VlnPlot(seu, features = "percent.mt",   pt.size = 0, ncol = 1) +
    ggtitle("Mitochondrial %") + NoLegend()

  p_violin <- p1 | p2 | p3
  ggsave(file.path(plot_dir, "01_qc_violins.pdf"), p_violin, width = 12, height = 5)

  p_scatter <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                               group.by = "sample") +
    ggtitle("UMI vs Genes") + NoLegend()
  ggsave(file.path(plot_dir, "02_scatter_qc.pdf"), p_scatter, width = 6, height = 5)
}


# ==============================================================================
#  STEP 3 — DOUBLET DETECTION (scDblFinder)
# ==============================================================================

detect_doublets <- function(seu) {
  if (!RUN_SCDBLFINDER) {
    seu$doublet_score     <- 0
    seu$predicted_doublet <- FALSE
    return(seu)
  }

  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce, dbr = EXPECTED_DBR, verbose = FALSE)

  seu$doublet_score     <- sce$scDblFinder.score
  seu$predicted_doublet <- sce$scDblFinder.class == "doublet"

  n_dbl <- sum(seu$predicted_doublet)
  message(sprintf("Doublets detected: %d (%.1f%%)", n_dbl, 100 * n_dbl / ncol(seu)))
  seu
}


# ==============================================================================
#  STEP 4 — FILTER CELLS AND GENES
# ==============================================================================

filter_cells_and_genes <- function(seu) {
  n_before <- ncol(seu)

  keep_cells <- (
    seu$nFeature_RNA >= MIN_GENES &
    seu$nFeature_RNA <= MAX_GENES &
    seu$nCount_RNA   >= MIN_COUNTS &
    seu$nCount_RNA   <= MAX_COUNTS &
    seu$percent.mt   <= MAX_MT_PCT
  )
  if (RUN_SCDBLFINDER) keep_cells <- keep_cells & !seu$predicted_doublet

  seu <- seu[, keep_cells]

  # Gene filter
  gene_counts <- rowSums(GetAssayData(seu, layer = "counts") > 0)
  keep_genes  <- gene_counts >= MIN_CELLS_PER_GENE
  seu         <- seu[keep_genes, ]

  message(sprintf(
    "Cells:  %d -> %d (%d removed)",
    n_before, ncol(seu), n_before - ncol(seu)
  ))
  message(sprintf("Genes expressed in >= %d cells retained: %d",
                  MIN_CELLS_PER_GENE, nrow(seu)))
  seu
}

save_before_after_plots <- function(seu_raw, seu_filtered, plot_dir) {
  make_df <- function(s, label) {
    data.frame(
      nCount_RNA   = s$nCount_RNA,
      nFeature_RNA = s$nFeature_RNA,
      percent.mt   = s$percent.mt,
      group        = label
    )
  }
  df <- rbind(make_df(seu_raw, "Before"), make_df(seu_filtered, "After"))

  p1 <- ggplot(df, aes(group, nCount_RNA, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    scale_fill_manual(values = c("Before" = "#AEC6E8", "After" = "#A8D5A2")) +
    theme_classic() + NoLegend() + ggtitle("UMI Counts")

  p2 <- ggplot(df, aes(group, nFeature_RNA, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    scale_fill_manual(values = c("Before" = "#AEC6E8", "After" = "#A8D5A2")) +
    theme_classic() + NoLegend() + ggtitle("Genes Detected")

  p3 <- ggplot(df, aes(group, percent.mt, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    scale_fill_manual(values = c("Before" = "#AEC6E8", "After" = "#A8D5A2")) +
    theme_classic() + NoLegend() + ggtitle("MT %")

  ggsave(file.path(plot_dir, "03_before_after_filter.pdf"),
         p1 | p2 | p3, width = 12, height = 5)
}


# ==============================================================================
#  MAIN
# ==============================================================================

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(OUTPUT_DIR, "plots")

  seu_raw      <- load_data()
  seu_metrics  <- compute_qc_metrics(seu_raw)
  save_qc_plots(seu_metrics, plot_dir)
  seu_scrubbed <- detect_doublets(seu_metrics)
  seu_filtered <- filter_cells_and_genes(seu_scrubbed)
  save_before_after_plots(seu_metrics, seu_filtered, plot_dir)

  out_path <- file.path(OUTPUT_DIR, paste0(SAMPLE_NAME, "_qc.rds"))
  saveRDS(seu_filtered, out_path)
  message("Saved -> ", out_path)
}

main()
