# 03_cell_annotation.R
# =====================
# Assign cell type identities to clusters.
#
# WHAT THIS SCRIPT DOES:
#   1. SingleR — automated reference-based annotation
#   2. Marker gene scoring — score each cell against canonical gene sets
#   3. Cluster marker discovery — Wilcoxon test to find top DE genes per cluster
#   4. Print a per-cluster consensus table so you can review and finalise labels
#   5. Apply the final annotation to seurat$cell_type
#   6. Save the annotated Seurat object
#
# HOW TO USE:
#   1. Edit CONFIGURATION below.
#   2. Run:  Rscript 03_cell_annotation.R
#   3. Review the CLUSTER ANNOTATION TABLE printed to the console.
#   4. Edit MANUAL_ANNOTATION to override any automated label.
#   5. Re-run to apply your manual labels and save.
#
# OUTPUT:
#   results/03_annotation/seurat_annotated.rds
#   results/03_annotation/cluster_markers.csv
#   results/03_annotation/plots/

suppressMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH <- "results/02_integration/seurat_integrated.rds"

# The cluster column to annotate (produced by script 02).
CLUSTER_KEY <- "leiden_0.50"

# "human" or "mouse"
ORGANISM <- "human"

# --- SingleR reference --------------------------------------------------------
# Options (celldex):
#   "HumanPrimaryCellAtlasData"   broad human cell types
#   "BlueprintEncodeData"         immune and stromal human cells
#   "ImmGenData"                  mouse immune cells
#   "MouseRNAseqData"             broad mouse cell types
SINGLER_REFERENCE <- "HumanPrimaryCellAtlasData"

# Use fine-grained labels ("label.fine") or broad labels ("label.main").
SINGLER_LABEL_TYPE <- "label.main"

# --- Cluster markers ----------------------------------------------------------
N_TOP_MARKER_GENES <- 10L

# --- Manual annotation override -----------------------------------------------
# After reviewing the printed table, add entries here to override automated labels.
# Format:  "cluster_id" = "Your Cell Type Name"
# Leave as an empty list to use SingleR labels for all clusters.
MANUAL_ANNOTATION <- list()
# Example:
# MANUAL_ANNOTATION <- list("0" = "CD4+ T cell", "3" = "Monocyte")

# --- Output -------------------------------------------------------------------
OUTPUT_DIR <- "results/03_annotation"


# ==============================================================================
#  HELPERS — canonical marker genes
# ==============================================================================

get_marker_genes <- function(organism) {
  if (tolower(organism) == "mouse") {
    list(
      "T cell"          = c("Cd3d", "Cd3e", "Cd3g"),
      "CD4+ T cell"     = c("Cd4", "Cd40lg", "Foxp3"),
      "CD8+ T cell"     = c("Cd8a", "Cd8b1", "Gzmb"),
      "NK cell"         = c("Nkg7", "Klrb1c", "Ncr1"),
      "B cell"          = c("Cd19", "Ms4a1", "Cd79a"),
      "Plasma cell"     = c("Sdc1", "Mzb1", "Xbp1"),
      "Monocyte"        = c("Cd14", "Csf1r", "Fcgr3"),
      "Macrophage"      = c("C1qa", "C1qb", "Adgre1"),
      "Dendritic cell"  = c("Itgax", "Siglech", "Ccr7"),
      "Neutrophil"      = c("S100a8", "S100a9", "Csf3r"),
      "Mast cell"       = c("Cpa3", "Mcpt4", "Fcer1a"),
      "Epithelial"      = c("Epcam", "Krt8", "Krt18"),
      "Fibroblast"      = c("Col1a1", "Col3a1", "Pdgfra"),
      "Endothelial"     = c("Pecam1", "Cdh5", "Kdr"),
      "Erythrocyte"     = c("Hba-a1", "Hba-a2", "Gypa")
    )
  } else {
    list(
      "T cell"          = c("CD3D", "CD3E", "CD3G"),
      "CD4+ T cell"     = c("CD4", "CD40LG", "FOXP3"),
      "CD8+ T cell"     = c("CD8A", "CD8B", "GZMB"),
      "NK cell"         = c("NKG7", "KLRB1", "NCR1"),
      "B cell"          = c("CD19", "MS4A1", "CD79A"),
      "Plasma cell"     = c("SDC1", "MZB1", "XBP1"),
      "Monocyte"        = c("CD14", "CSF1R", "FCGR3A"),
      "Macrophage"      = c("C1QA", "C1QB", "MRC1"),
      "Dendritic cell"  = c("ITGAX", "SIGLEC6", "CCR7"),
      "Neutrophil"      = c("S100A8", "S100A9", "CSF3R"),
      "Mast cell"       = c("CPA3", "TPSAB1", "FCER1A"),
      "Epithelial"      = c("EPCAM", "KRT8", "KRT18"),
      "Fibroblast"      = c("COL1A1", "COL3A1", "PDGFRA"),
      "Endothelial"     = c("PECAM1", "CDH5", "KDR"),
      "Erythrocyte"     = c("HBA1", "HBA2", "GYPA")
    )
  }
}


# ==============================================================================
#  STEP 1 — LOAD DATA
# ==============================================================================

load_data <- function() {
  if (!file.exists(INPUT_PATH)) stop("File not found: ", INPUT_PATH)
  seu <- readRDS(INPUT_PATH)

  if (!CLUSTER_KEY %in% colnames(seu@meta.data)) {
    available <- grep("leiden_", colnames(seu@meta.data), value = TRUE)
    stop(sprintf("Column '%s' not found.\nAvailable: %s",
                 CLUSTER_KEY, paste(available, collapse = ", ")))
  }
  seu
}


# ==============================================================================
#  STEP 2 — SINGLER ANNOTATION
# ==============================================================================

run_singler <- function(seu) {
  ref_fn <- tryCatch(
    get(SINGLER_REFERENCE, envir = asNamespace("celldex")),
    error = function(e) {
      stop("Unknown SingleR reference: ", SINGLER_REFERENCE,
           "\nOptions: HumanPrimaryCellAtlasData, BlueprintEncodeData, ",
           "ImmGenData, MouseRNAseqData")
    }
  )
  ref <- ref_fn()

  counts <- GetAssayData(seu, layer = "data")
  pred <- SingleR(
    test      = counts,
    ref       = ref,
    labels    = ref[[SINGLER_LABEL_TYPE]],
    de.method = "wilcox"
  )

  seu$singler_label <- pred$labels
  seu$singler_score <- apply(pred$scores, 1, max)
  message(sprintf("SingleR annotation complete. Top labels: %s",
                  paste(head(sort(table(pred$labels), decreasing = TRUE), 5L), collapse = " | ")))
  seu
}


# ==============================================================================
#  STEP 3 — MARKER GENE SCORING
# ==============================================================================

run_marker_scoring <- function(seu) {
  markers <- get_marker_genes(ORGANISM)
  genes_in_data <- rownames(seu)

  scored_cts <- character(0)
  for (ct in names(markers)) {
    present <- intersect(markers[[ct]], genes_in_data)
    if (length(present) < 2L) next
    score_key <- gsub("[ +/]", "_", ct)
    seu <- AddModuleScore(seu, features = list(present),
                          name = paste0("score_", score_key), seed = 42)
    scored_cts <- c(scored_cts, ct)
  }

  score_cols <- grep("^score_", colnames(seu@meta.data), value = TRUE)
  if (length(score_cols) > 0L) {
    score_mat <- as.matrix(seu@meta.data[, score_cols, drop = FALSE])
    best_idx  <- apply(score_mat, 1, which.max)
    seu$marker_score_label <- scored_cts[best_idx]
  } else {
    seu$marker_score_label <- "unknown"
    message("Warning: no marker genes overlapped with dataset. Check ORGANISM setting.")
  }
  seu
}


# ==============================================================================
#  STEP 4 — DATA-DRIVEN CLUSTER MARKERS
# ==============================================================================

find_cluster_markers <- function(seu, out_dir) {
  Idents(seu) <- CLUSTER_KEY
  markers <- FindAllMarkers(
    seu,
    only.pos  = TRUE,
    min.pct   = 0.25,
    logfc.threshold = 0.25,
    test.use  = "wilcox",
    verbose   = FALSE
  )
  markers <- markers |>
    group_by(cluster) |>
    slice_max(order_by = avg_log2FC, n = N_TOP_MARKER_GENES) |>
    ungroup()

  write.csv(markers, file.path(out_dir, "cluster_markers.csv"), row.names = FALSE)
  seu
}


# ==============================================================================
#  STEP 5 — CONSENSUS TABLE & FINAL ANNOTATION
# ==============================================================================

print_annotation_table <- function(seu) {
  meta <- seu@meta.data
  clusters <- sort(unique(as.character(meta[[CLUSTER_KEY]])),
                   method = "radix")

  cat(sprintf("\n%-10s %7s  %-30s  %-30s\n",
              "Cluster", "Cells", "SingleR", "Marker Score"))
  cat(strrep("-", 82), "\n")

  for (cl in clusters) {
    mask <- as.character(meta[[CLUSTER_KEY]]) == cl
    sr_label <- names(sort(table(meta$singler_label[mask]),
                           decreasing = TRUE))[1]
    ms_label <- names(sort(table(meta$marker_score_label[mask]),
                           decreasing = TRUE))[1]
    n <- sum(mask)
    override <- if (!is.null(MANUAL_ANNOTATION[[cl]])) {
      sprintf("  <- MANUAL: %s", MANUAL_ANNOTATION[[cl]])
    } else ""
    cat(sprintf("%-10s %7d  %-30s  %-30s%s\n",
                cl, n, sr_label, ms_label, override))
  }
  cat("\nTo override a label, add it to MANUAL_ANNOTATION and re-run.\n\n")
}

apply_annotation <- function(seu) {
  meta      <- seu@meta.data
  clusters  <- unique(as.character(meta[[CLUSTER_KEY]]))
  mapping   <- setNames(character(length(clusters)), clusters)

  for (cl in clusters) {
    if (!is.null(MANUAL_ANNOTATION[[cl]])) {
      mapping[cl] <- MANUAL_ANNOTATION[[cl]]
    } else if ("singler_label" %in% colnames(meta)) {
      mask <- as.character(meta[[CLUSTER_KEY]]) == cl
      mapping[cl] <- names(sort(table(meta$singler_label[mask]),
                                decreasing = TRUE))[1]
    } else {
      mapping[cl] <- paste0("Cluster_", cl)
    }
  }

  seu$cell_type <- factor(mapping[as.character(meta[[CLUSTER_KEY]])])
  seu
}

save_umap_plot <- function(seu, plot_dir) {
  p <- DimPlot(seu, reduction = "umap", group.by = "cell_type",
               label = TRUE, label.size = 3.5, pt.size = 0.3, repel = TRUE) +
    ggtitle("UMAP — Final Cell Type Annotation") +
    theme_classic() +
    NoLegend()
  ggsave(file.path(plot_dir, "01_umap_cell_types.pdf"), p, width = 9, height = 8)
}


# ==============================================================================
#  MAIN
# ==============================================================================

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(OUTPUT_DIR, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  seu <- load_data()
  seu <- run_singler(seu)
  seu <- run_marker_scoring(seu)
  seu <- find_cluster_markers(seu, OUTPUT_DIR)
  print_annotation_table(seu)
  seu <- apply_annotation(seu)
  save_umap_plot(seu, plot_dir)

  out_path <- file.path(OUTPUT_DIR, "seurat_annotated.rds")
  saveRDS(seu, out_path)
  message("Saved -> ", out_path)
}

main()
