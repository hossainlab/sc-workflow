# 05_trajectory_analysis.R
# ==========================
# Trajectory and pseudotime analysis using Monocle 3.
#
# WHAT THIS SCRIPT DOES:
#   1. Load annotated data and subset to selected cell types / lineage
#   2. Convert Seurat object to Monocle 3 cell_data_set (CDS)
#   3. Learn the principal graph (trajectory)
#   4. Order cells in pseudotime
#   5. Identify genes that change along pseudotime (graph_test)
#   6. Save trajectory CDS and pseudotime-annotated Seurat object
#
# WHEN TO USE:
#   Developmental systems (stem cell differentiation), time-course experiments
#   (treatment response), or any biology with cell state transitions.
#
# HOW TO USE:
#   1. Edit CONFIGURATION below.
#   2. Run:  Rscript 05_trajectory_analysis.R
#
# OUTPUT:
#   results/05_trajectory/seurat_trajectory.rds
#   results/05_trajectory/pseudotime_genes.csv
#   results/05_trajectory/plots/

suppressMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH    <- "results/03_annotation/seurat_annotated.rds"
CELL_TYPE_COL <- "cell_type"

# --- Lineage selection --------------------------------------------------------
# Cell types to include in the trajectory.
# Leave as character(0) to use all cell types.
SUBSET_CELL_TYPES <- character(0)
# Example: c("Stem cell", "Progenitor", "Effector T cell")

# Cell type to use as the root (start) of the trajectory.
# Leave as "" for automatic detection.
ROOT_CELL_TYPE <- ""

# --- Monocle 3 parameters -----------------------------------------------------
# Number of partitions for learn_graph (set FALSE to use all cells as one partition).
USE_PARTITION <- FALSE

# Fraction of cells to use per cell type (reduce for large datasets).
DOWNSAMPLE_FRAC <- 1.0   # set to e.g. 0.3 to use 30% of cells

# Minimum cells per cell type after downsampling.
MIN_CELLS_KEPT <- 50L

# Number of top pseudotime-variable genes to report.
N_TOP_PT_GENES <- 200L

# --- Output -------------------------------------------------------------------
OUTPUT_DIR  <- "results/05_trajectory"
OUTPUT_PATH <- file.path(OUTPUT_DIR, "seurat_trajectory.rds")


# ==============================================================================
#  STEP 1 — LOAD & SUBSET
# ==============================================================================

load_and_subset <- function() {
  if (!file.exists(INPUT_PATH)) stop("File not found: ", INPUT_PATH)
  seu <- readRDS(INPUT_PATH)

  if (length(SUBSET_CELL_TYPES) > 0L) {
    mask <- seu@meta.data[[CELL_TYPE_COL]] %in% SUBSET_CELL_TYPES
    seu  <- seu[, mask]
    message(sprintf("Subsetted to %d cell types: %d cells",
                    length(SUBSET_CELL_TYPES), ncol(seu)))
  }

  if (DOWNSAMPLE_FRAC < 1.0) {
    ct_vec    <- as.character(seu@meta.data[[CELL_TYPE_COL]])
    keep_cells <- unlist(tapply(seq_len(ncol(seu)), ct_vec, function(idx) {
      n_keep <- max(MIN_CELLS_KEPT, floor(length(idx) * DOWNSAMPLE_FRAC))
      n_keep <- min(n_keep, length(idx))
      sample(idx, n_keep)
    }))
    seu <- seu[, sort(keep_cells)]
    message(sprintf("Downsampled to %d cells (frac = %.2f)", ncol(seu), DOWNSAMPLE_FRAC))
  }

  seu
}


# ==============================================================================
#  STEP 2 — SEURAT -> MONOCLE 3 CONVERSION
# ==============================================================================

seurat_to_cds <- function(seu) {
  # Monocle 3 requires raw counts, gene metadata, and cell metadata.
  counts    <- GetAssayData(seu, layer = "counts")
  cell_meta <- seu@meta.data

  gene_meta <- data.frame(
    gene_short_name = rownames(counts),
    row.names       = rownames(counts)
  )

  cds <- new_cell_data_set(
    expression_data = counts,
    cell_metadata   = cell_meta,
    gene_metadata   = gene_meta
  )

  # Transfer UMAP embedding from Seurat so we can reuse it.
  if ("umap" %in% names(seu@reductions)) {
    reducedDims(cds)[["UMAP"]] <- Embeddings(seu, "umap")
  }

  cds
}


# ==============================================================================
#  STEP 3 — PREPROCESS, LEARN GRAPH & ORDER CELLS
# ==============================================================================

run_monocle3 <- function(cds, seu) {
  # Preprocessing (uses PCA stored in Seurat if available, otherwise recomputes)
  if ("pca" %in% names(seu@reductions)) {
    reducedDims(cds)[["PCA"]] <- Embeddings(seu, "pca")
  } else {
    cds <- preprocess_cds(cds, num_dim = 50L, verbose = FALSE)
  }

  # UMAP (reuse Seurat UMAP so cell positions are consistent across scripts)
  if (!"UMAP" %in% names(reducedDims(cds))) {
    cds <- reduce_dimension(cds, reduction_method = "UMAP", verbose = FALSE)
  }

  # Cluster cells (required internally by learn_graph)
  cds <- cluster_cells(cds, verbose = FALSE)

  # Learn the principal graph
  cds <- learn_graph(cds, use_partition = USE_PARTITION, verbose = FALSE)

  # Set root
  if (nchar(ROOT_CELL_TYPE) > 0L) {
    root_mask <- colData(cds)[[CELL_TYPE_COL]] == ROOT_CELL_TYPE
    if (!any(root_mask)) {
      message("ROOT_CELL_TYPE not found — using automatic root detection.")
      cds <- order_cells(cds, verbose = FALSE)
    } else {
      root_cell <- colnames(cds)[root_mask][1]
      cds <- order_cells(cds, root_cells = root_cell, verbose = FALSE)
    }
  } else {
    cds <- order_cells(cds, verbose = FALSE)
  }

  # Copy pseudotime back to Seurat
  seu$pseudotime <- pseudotime(cds, reduction_method = "UMAP")
  list(cds = cds, seu = seu)
}


# ==============================================================================
#  STEP 4 — PSEUDOTIME-VARIABLE GENES
# ==============================================================================

find_pseudotime_genes <- function(cds) {
  message("Running graph_test to find pseudotime-variable genes ...")
  pt_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1L)
  pt_genes <- pt_genes[order(pt_genes$q_value), ]
  head(pt_genes, N_TOP_PT_GENES)
}


# ==============================================================================
#  STEP 5 — VISUALISATION
# ==============================================================================

save_trajectory_plots <- function(seu, cds, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Pseudotime on UMAP (from Seurat)
  p1 <- FeaturePlot(seu, features = "pseudotime", reduction = "umap") +
    scale_colour_viridis_c() +
    ggtitle("UMAP — Diffusion Pseudotime") +
    theme_classic()
  ggsave(file.path(plot_dir, "01_umap_pseudotime.pdf"), p1, width = 7, height = 6)

  # Cell types on UMAP
  p2 <- DimPlot(seu, reduction = "umap", group.by = CELL_TYPE_COL,
                label = TRUE, pt.size = 0.4, repel = TRUE) +
    ggtitle("UMAP — Cell Types") +
    theme_classic() + NoLegend()
  ggsave(file.path(plot_dir, "02_umap_celltypes.pdf"), p2, width = 7, height = 6)

  # Monocle 3 trajectory plot
  p3 <- plot_cells(
    cds,
    color_cells_by   = CELL_TYPE_COL,
    label_cell_groups = FALSE,
    label_leaves     = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 3
  ) + ggtitle("Monocle 3 Trajectory")
  ggsave(file.path(plot_dir, "03_monocle3_trajectory.pdf"), p3, width = 8, height = 7)
}


# ==============================================================================
#  MAIN
# ==============================================================================

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(OUTPUT_DIR, "plots")

  seu          <- load_and_subset()
  cds          <- seurat_to_cds(seu)
  result       <- run_monocle3(cds, seu)
  cds          <- result$cds
  seu          <- result$seu
  pt_genes     <- find_pseudotime_genes(cds)
  save_trajectory_plots(seu, cds, plot_dir)

  write.csv(pt_genes, file.path(OUTPUT_DIR, "pseudotime_genes.csv"))
  saveRDS(seu, OUTPUT_PATH)
  message("Saved -> ", OUTPUT_PATH)
}

main()
