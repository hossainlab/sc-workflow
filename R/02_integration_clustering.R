# 02_integration_clustering.R
# ============================
# Merge QC-filtered samples, correct batch effects, and cluster cells.
#
# WHAT THIS SCRIPT DOES:
#   1. Load one or more QC-filtered .rds files from script 01
#   2. Normalize counts and select highly variable genes (per sample)
#   3. Merge samples and run PCA
#   4. Apply Harmony batch correction
#   5. Build k-nearest-neighbour graph and embed in UMAP
#   6. Leiden/Louvain clustering at multiple resolutions
#   7. Save integrated Seurat object
#
# BEST PRACTICE:
#   Normalize and find HVGs independently per sample before merging.
#   Never normalize the merged object globally first.
#
# HOW TO USE:
#   1. Edit CONFIGURATION below — especially INPUT_PATHS and BATCH_KEY.
#   2. Run:  Rscript 02_integration_clustering.R
#
# OUTPUT:
#   results/02_integration/seurat_integrated.rds
#   results/02_integration/plots/

suppressMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

# Named list or character vector of .rds files from script 01.
INPUT_PATHS <- c(
  "results/01_qc/sample1_qc.rds"
  # "results/01_qc/sample2_qc.rds"
)

# Column in Seurat metadata that identifies which sample each cell came from.
BATCH_KEY <- "sample"

# --- Normalisation & HVG parameters ------------------------------------------
# Normalise each cell to this total count before log-transformation.
TARGET_SUM <- 1e4

# Number of highly variable genes to select per sample.
N_HVGS <- 3000L

# Number of principal components to compute and use.
N_PCS <- 50L

# --- Harmony batch correction -------------------------------------------------
RUN_HARMONY   <- TRUE
HARMONY_THETA <- 2.0     # diversity penalty (higher = stronger correction)
HARMONY_ITERS <- 10L     # number of Harmony iterations

# --- Neighbourhood graph & clustering -----------------------------------------
N_NEIGHBOURS <- 30L
UMAP_MIN_DIST <- 0.3

# Clustering resolutions (higher = more, smaller clusters).
CLUSTERING_RESOLUTIONS <- c(0.3, 0.5, 0.8, 1.0)

# Resolution used for display plots (must be in CLUSTERING_RESOLUTIONS).
DEFAULT_RESOLUTION <- 0.5

# --- Output -------------------------------------------------------------------
OUTPUT_DIR <- "results/02_integration"


# ==============================================================================
#  STEP 1 — LOAD & MERGE SAMPLES
# ==============================================================================

load_and_merge <- function() {
  seu_list <- lapply(INPUT_PATHS, function(p) {
    if (!file.exists(p)) stop("File not found: ", p)
    s <- readRDS(p)
    if (!BATCH_KEY %in% colnames(s@meta.data)) {
      s[[BATCH_KEY]] <- tools::file_path_sans_ext(basename(p))
    }
    s
  })

  if (length(seu_list) == 1L) return(seu_list[[1]])

  merged <- merge(
    x    = seu_list[[1]],
    y    = seu_list[-1],
    add.cell.ids = sapply(seu_list, function(s) unique(s[[BATCH_KEY]][[1]])),
    merge.data   = FALSE
  )
  merged
}


# ==============================================================================
#  STEP 2 — NORMALISATION, HVG SELECTION & PCA
# ==============================================================================

normalise_and_pca <- function(seu) {
  # Store raw counts for DE analysis in script 04
  seu <- JoinLayers(seu)
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu[[BATCH_KEY]])

  samples <- unique(seu[[BATCH_KEY]][[1]])

  if (length(samples) > 1L) {
    # Normalise per sample, find HVGs per sample, take union of shared HVGs.
    # HVG is "shared" if selected in >= 2 samples (or all samples if 2 total).
    hvg_list <- lapply(samples, function(s) {
      sub <- subset(seu, cells = which(seu[[BATCH_KEY]][[1]] == s))
      sub <- NormalizeData(sub, normalization.method = "LogNormalize",
                           scale.factor = TARGET_SUM, verbose = FALSE)
      sub <- FindVariableFeatures(sub, selection.method = "vst",
                                  nfeatures = N_HVGS, verbose = FALSE)
      VariableFeatures(sub)
    })

    hvg_table     <- table(unlist(hvg_list))
    min_samples   <- min(2L, length(samples))
    shared_hvgs   <- names(hvg_table[hvg_table >= min_samples])

    seu <- JoinLayers(seu)
    seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                         scale.factor = TARGET_SUM, verbose = FALSE)
    VariableFeatures(seu) <- shared_hvgs
    message(sprintf("Shared HVGs across >= %d samples: %d", min_samples, length(shared_hvgs)))

  } else {
    seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                         scale.factor = TARGET_SUM, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst",
                                nfeatures = N_HVGS, verbose = FALSE)
  }

  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = N_PCS, seed.use = 42, verbose = FALSE)
  seu
}

save_elbow_plot <- function(seu, plot_dir) {
  p <- ElbowPlot(seu, ndims = N_PCS) +
    ggtitle("PCA Elbow Plot") +
    theme_classic()
  ggsave(file.path(plot_dir, "01_pca_elbow.pdf"), p, width = 7, height = 4)
}


# ==============================================================================
#  STEP 3 — HARMONY BATCH CORRECTION
# ==============================================================================

run_harmony <- function(seu) {
  if (!RUN_HARMONY) return(list(seu = seu, reduction = "pca"))

  seu <- RunHarmony(
    seu,
    group.by.vars   = BATCH_KEY,
    theta           = HARMONY_THETA,
    max.iter.harmony = HARMONY_ITERS,
    reduction       = "pca",
    reduction.save  = "harmony",
    verbose         = FALSE
  )
  list(seu = seu, reduction = "harmony")
}


# ==============================================================================
#  STEP 4 — UMAP & CLUSTERING
# ==============================================================================

embed_and_cluster <- function(seu, reduction) {
  seu <- FindNeighbors(seu, reduction = reduction, dims = 1:N_PCS,
                       k.param = N_NEIGHBOURS, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = reduction, dims = 1:N_PCS,
                 min.dist = UMAP_MIN_DIST, seed.use = 42, verbose = FALSE)

  for (res in CLUSTERING_RESOLUTIONS) {
    col <- sprintf("leiden_%.2f", res)
    seu <- FindClusters(seu, resolution = res, algorithm = 4L,
                        random.seed = 42, cluster.name = col, verbose = FALSE)
  }
  seu
}

save_umap_plots <- function(seu, plot_dir) {
  col <- sprintf("leiden_%.2f", DEFAULT_RESOLUTION)
  if (!col %in% colnames(seu@meta.data)) {
    col <- grep("leiden_", colnames(seu@meta.data), value = TRUE)[1]
  }

  p_cluster <- DimPlot(seu, reduction = "umap", group.by = col,
                        label = TRUE, pt.size = 0.3) +
    ggtitle(sprintf("UMAP — Leiden clusters (res %.1f)", DEFAULT_RESOLUTION)) +
    theme_classic()

  p_batch <- DimPlot(seu, reduction = "umap", group.by = BATCH_KEY,
                      pt.size = 0.3) +
    ggtitle(sprintf("UMAP — Batch (%s)", BATCH_KEY)) +
    theme_classic()

  ggsave(file.path(plot_dir, "02_umap_clusters.pdf"), p_cluster, width = 8, height = 7)
  ggsave(file.path(plot_dir, "03_umap_batch.pdf"),    p_batch,   width = 8, height = 7)
}


# ==============================================================================
#  MAIN
# ==============================================================================

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(OUTPUT_DIR, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  seu_merged           <- load_and_merge()
  seu_pca              <- normalise_and_pca(seu_merged)
  save_elbow_plot(seu_pca, plot_dir)
  res                  <- run_harmony(seu_pca)
  seu_integrated       <- embed_and_cluster(res$seu, res$reduction)
  save_umap_plots(seu_integrated, plot_dir)

  out_path <- file.path(OUTPUT_DIR, "seurat_integrated.rds")
  saveRDS(seu_integrated, out_path)
  message("Saved -> ", out_path)
}

main()
