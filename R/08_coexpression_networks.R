# 08_coexpression_networks.R
# ============================
# Gene co-expression network analysis using WGCNA.
#
# WHAT THIS SCRIPT DOES:
#   1. Subset to target cell type and build metacells (k-NN aggregation)
#      — metacells average k nearest neighbours to reduce technical zeros
#   2. Soft-thresholding power selection (scale-free topology test)
#   3. WGCNA network construction and module detection
#   4. Hub gene identification per module (by intramodular connectivity kME)
#   5. Pathway enrichment per module via enrichR
#   6. Save all results
#
# WHY METACELLS?
#   Single-cell data is 80-90% sparse. Averaging k nearest neighbours
#   in PCA space fills technical zeros without losing cell-type specificity.
#
# HOW TO USE:
#   1. Edit CONFIGURATION below.
#   2. Run:  Rscript 08_coexpression_networks.R
#
# OUTPUT:
#   results/08_coexpression/module_assignments.csv
#   results/08_coexpression/module_eigengenes.csv
#   results/08_coexpression/hub_genes.csv
#   results/08_coexpression/module_enrichment.csv
#   results/08_coexpression/plots/

suppressMessages({
  library(Seurat)
  library(WGCNA)
  library(FNN)
  library(ggplot2)
  library(dplyr)
  library(enrichR)
})

options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads()


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH    <- "results/03_annotation/seurat_annotated.rds"
CELL_TYPE_COL <- "cell_type"
SAMPLE_COL    <- "sample"
CONDITION_COL <- "condition"    # used for module-trait correlation

# --- Target cell type ---------------------------------------------------------
# Leave "" to use all cells (not recommended — very slow and noisy).
TARGET_CELL_TYPE <- "CD4+ T cell"

# --- Metacell parameters ------------------------------------------------------
METACELL_K        <- 20L     # neighbours to average per metacell
MAX_SHARED        <- 10L     # max cell overlap between two metacells
MIN_EXPR_FRACTION <- 0.05    # keep genes expressed in >= 5% of metacells

# --- WGCNA parameters ---------------------------------------------------------
# Soft-thresholding power β.
# Set to 0 to auto-detect (chooses lowest β where scale-free R² >= 0.85).
SOFT_POWER     <- 0L
MIN_MODULE_SIZE <- 30L
MERGE_CUT_HEIGHT <- 0.25   # modules with correlation > (1 - 0.25) are merged

# --- Pathway enrichment -------------------------------------------------------
RUN_ENRICHR        <- TRUE
MODULE_TOP_N_GENES <- 50L    # top genes per module submitted to enrichR
ENRICHR_LIBS       <- c("MSigDB_Hallmark_2020", "KEGG_2021_Human")

# --- Output -------------------------------------------------------------------
OUTPUT_DIR <- "results/08_coexpression"


# ==============================================================================
#  HELPERS
# ==============================================================================

header <- function(title) {
  message()
  message(strrep("=", 60))
  message("  ", title)
  message(strrep("=", 60))
}


# ==============================================================================
#  STEP 1 — LOAD & SUBSET
# ==============================================================================

load_and_subset <- function() {
  header("Step 1 — Load & Subset Data")
  if (!file.exists(INPUT_PATH)) stop("File not found: ", INPUT_PATH)
  seu <- readRDS(INPUT_PATH)
  message(sprintf("  Loaded %d cells x %d genes", ncol(seu), nrow(seu)))

  if (nchar(TARGET_CELL_TYPE) > 0L) {
    mask <- seu@meta.data[[CELL_TYPE_COL]] == TARGET_CELL_TYPE
    seu  <- seu[, mask]
    if (ncol(seu) < 100L) {
      stop(sprintf("Only %d cells for '%s'. Need >= 100 cells.",
                   ncol(seu), TARGET_CELL_TYPE))
    }
    message(sprintf("  Subset to '%s': %d cells", TARGET_CELL_TYPE, ncol(seu)))
  } else {
    message("  Using all cells (TARGET_CELL_TYPE is empty)")
  }
  seu
}


# ==============================================================================
#  STEP 2 — METACELL CONSTRUCTION
# ==============================================================================

build_metacells <- function(seu) {
  header("Step 2 — Metacell Construction")
  message(sprintf("  k = %d nearest neighbours per metacell", METACELL_K))

  # Use raw counts
  counts_mat <- as.matrix(GetAssayData(seu, layer = "counts"))   # genes x cells
  counts_mat <- t(counts_mat)                                     # cells x genes

  # PCA embedding for kNN
  if ("pca" %in% names(seu@reductions)) {
    pca_mat <- Embeddings(seu, "pca")
  } else {
    pca_mat <- counts_mat[, seq_len(min(50L, ncol(counts_mat)))]
  }

  samples <- if (SAMPLE_COL %in% colnames(seu@meta.data)) {
    unique(as.character(seu@meta.data[[SAMPLE_COL]]))
  } else {
    "all"
  }

  mc_rows <- list()
  mc_meta <- list()

  for (s in samples) {
    if (s == "all") {
      idx <- seq_len(nrow(pca_mat))
    } else {
      idx <- which(seu@meta.data[[SAMPLE_COL]] == s)
    }

    pca_s    <- pca_mat[idx, , drop = FALSE]
    counts_s <- counts_mat[idx, , drop = FALSE]

    # kNN
    nn <- FNN::get.knn(pca_s, k = METACELL_K)$nn.index

    used      <- logical(nrow(pca_s))
    n_target  <- max(1L, nrow(pca_s) %/% METACELL_K)

    for (i in seq_len(nrow(pca_s))) {
      knn_idx <- nn[i, ]
      overlap <- sum(used[knn_idx])
      if (overlap > MAX_SHARED) next

      mc_rows <- c(mc_rows, list(colMeans(counts_s[knn_idx, , drop = FALSE])))
      mc_meta <- c(mc_meta, list(data.frame(sample = s)))
      used[knn_idx] <- TRUE
      if (length(mc_rows) >= n_target) break
    }
  }

  if (length(mc_rows) == 0L) stop("Metacell construction failed. Try reducing METACELL_K.")

  mc_matrix <- do.call(rbind, mc_rows)   # metacells x genes
  mc_df     <- data.frame(do.call(rbind, mc_meta))

  # Gene filter: expressed in >= MIN_EXPR_FRACTION of metacells
  expr_frac <- colMeans(mc_matrix > 0)
  keep      <- expr_frac >= MIN_EXPR_FRACTION
  mc_matrix <- mc_matrix[, keep, drop = FALSE]

  message(sprintf("  Built %d metacells x %d genes", nrow(mc_matrix), ncol(mc_matrix)))
  list(matrix = mc_matrix, meta = mc_df)
}


# ==============================================================================
#  STEP 3 — SOFT-THRESHOLDING POWER SELECTION
# ==============================================================================

select_soft_power <- function(mc_matrix, plot_dir) {
  header("Step 3 — Soft-Thresholding Power Selection")

  if (SOFT_POWER > 0L) {
    message(sprintf("  Using manually set SOFT_POWER = %d", SOFT_POWER))
    return(SOFT_POWER)
  }

  powers  <- 1:30
  sft     <- pickSoftThreshold(mc_matrix, powerVector = powers,
                                networkType = "unsigned", verbose = 0)

  # Plot
  df <- data.frame(Power = sft$fitIndices$Power,
                   R2    = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
                   MeanK = sft$fitIndices$mean.k.)

  p1 <- ggplot(df, aes(Power, R2)) +
    geom_line() + geom_point() +
    geom_hline(yintercept = 0.85, linetype = "dashed", colour = "red") +
    labs(title = "Scale-Free Fit (R²)", x = "Soft Power (β)", y = "R²") +
    theme_classic()

  p2 <- ggplot(df, aes(Power, MeanK)) +
    geom_line() + geom_point() +
    labs(title = "Mean Connectivity", x = "Soft Power (β)", y = "Mean k") +
    theme_classic()

  ggsave(file.path(plot_dir, "00_soft_power_selection.pdf"),
         p1 + p2, width = 10, height = 4)

  auto_power <- sft$powerEstimate
  if (is.na(auto_power)) {
    auto_power <- 6L
    message("  R² >= 0.85 not reached — defaulting to β = 6. Consider reviewing the plot.")
  }
  message(sprintf("  Auto-selected soft power β = %d", auto_power))
  as.integer(auto_power)
}


# ==============================================================================
#  STEP 4 — WGCNA MODULE DETECTION
# ==============================================================================

run_wgcna <- function(mc_matrix, selected_power) {
  header("Step 4 — WGCNA Network Construction & Module Detection")

  # blockwiseModules handles large gene sets by processing in blocks
  net <- blockwiseModules(
    datExpr          = mc_matrix,
    power            = selected_power,
    networkType      = "unsigned",
    TOMType          = "unsigned",
    minModuleSize    = MIN_MODULE_SIZE,
    mergeCutHeight   = MERGE_CUT_HEIGHT,
    reassignThreshold = 0,
    numericLabels    = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs         = FALSE,
    verbose          = 0
  )

  n_modules <- length(unique(net$colors)) - 1L  # exclude grey (unassigned)
  message(sprintf("  %d co-expression modules detected", n_modules))
  net
}


# ==============================================================================
#  STEP 5 — HUB GENES (kME)
# ==============================================================================

get_hub_genes <- function(mc_matrix, net) {
  header("Step 5 — Hub Genes (Intramodular Connectivity kME)")

  kme_mat <- signedKME(mc_matrix, net$MEs, outputColumnName = "kME")

  modules <- unique(net$colors[net$colors != "grey"])
  rows    <- list()

  for (mod in modules) {
    col     <- paste0("kME", mod)
    if (!col %in% colnames(kme_mat)) next
    gene_mask <- net$colors == mod
    mod_kme   <- kme_mat[gene_mask, col, drop = TRUE]
    top_idx   <- order(mod_kme, decreasing = TRUE)[seq_len(min(10L, sum(gene_mask)))]
    rows[[mod]] <- data.frame(
      gene   = rownames(kme_mat)[gene_mask][top_idx],
      module = mod,
      kME    = mod_kme[top_idx]
    )
  }

  hub_df <- bind_rows(rows)
  message("  Top 3 hub genes per module:")
  for (mod in modules) {
    top3 <- head(hub_df$gene[hub_df$module == mod], 3L)
    message(sprintf("    %s: %s", mod, paste(top3, collapse = ", ")))
  }
  hub_df
}


# ==============================================================================
#  STEP 6 — MODULE EIGENGENE PLOT
# ==============================================================================

save_eigengene_plot <- function(net, plot_dir) {
  me_corr <- cor(net$MEs)

  me_df <- as.data.frame(as.table(me_corr))
  colnames(me_df) <- c("Module1", "Module2", "Correlation")

  p <- ggplot(me_df, aes(Module1, Module2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "#4C72B0", mid = "white", high = "#C44E52",
                         midpoint = 0, limits = c(-1, 1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7)) +
    labs(title = "Module Eigengene Correlation Matrix",
         x = NULL, y = NULL)

  ggsave(file.path(plot_dir, "01_eigengene_correlation.pdf"), p, width = 8, height = 7)
}


# ==============================================================================
#  STEP 7 — PATHWAY ENRICHMENT PER MODULE
# ==============================================================================

run_module_enrichment <- function(mc_matrix, net) {
  header("Step 7 — Pathway Enrichment per Module (enrichR)")

  if (!RUN_ENRICHR) {
    message("  Skipped (RUN_ENRICHR = FALSE)")
    return(NULL)
  }

  modules <- unique(net$colors[net$colors != "grey"])
  results <- list()

  for (mod in modules) {
    gene_mask <- net$colors == mod
    genes     <- colnames(mc_matrix)[gene_mask]
    genes     <- head(genes, MODULE_TOP_N_GENES)
    if (length(genes) < 5L) next

    tryCatch({
      enr <- enrichr(genes, ENRICHR_LIBS)
      for (lib in names(enr)) {
        r          <- head(enr[[lib]], 5L)
        r$module   <- mod
        r$n_genes  <- length(genes)
        r$library  <- lib
        results[[paste(mod, lib)]] <- r
      }
      message(sprintf("  Module %s: enrichment done", mod))
    }, error = function(e) {
      message(sprintf("  Module %s: enrichment failed — %s", mod, conditionMessage(e)))
    })
  }

  if (length(results) == 0L) return(NULL)
  bind_rows(results)
}


# ==============================================================================
#  MAIN
# ==============================================================================

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(OUTPUT_DIR, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  message()
  message(strrep("=", 60))
  message("  scRNA-seq Co-expression Networks (WGCNA)")
  message(strrep("=", 60))
  message(sprintf("  Input      : %s", INPUT_PATH))
  message(sprintf("  Cell type  : %s", if (nchar(TARGET_CELL_TYPE) > 0L) TARGET_CELL_TYPE else "all cells"))
  message(sprintf("  Metacell k : %d", METACELL_K))
  message(sprintf("  Output     : %s", OUTPUT_DIR))

  seu            <- load_and_subset()
  mc             <- build_metacells(seu)
  power          <- select_soft_power(mc$matrix, plot_dir)
  net            <- run_wgcna(mc$matrix, power)
  hub_genes      <- get_hub_genes(mc$matrix, net)
  save_eigengene_plot(net, plot_dir)
  enrichment_df  <- run_module_enrichment(mc$matrix, net)

  header("Save")

  # Module assignments
  assign_df <- data.frame(gene = colnames(mc$matrix), module = net$colors)
  write.csv(assign_df,  file.path(OUTPUT_DIR, "module_assignments.csv"),  row.names = FALSE)
  message("  Module assignments -> module_assignments.csv")

  # Module eigengenes
  write.csv(net$MEs,    file.path(OUTPUT_DIR, "module_eigengenes.csv"))
  message("  Module eigengenes  -> module_eigengenes.csv")

  # Hub genes
  write.csv(hub_genes,  file.path(OUTPUT_DIR, "hub_genes.csv"),           row.names = FALSE)
  message("  Hub genes          -> hub_genes.csv")

  # Pathway enrichment
  if (!is.null(enrichment_df)) {
    write.csv(enrichment_df, file.path(OUTPUT_DIR, "module_enrichment.csv"), row.names = FALSE)
    message("  Module enrichment  -> module_enrichment.csv")
  }

  message()
  message("  Pipeline complete. All results in: ", OUTPUT_DIR)
  message()
}

main()
