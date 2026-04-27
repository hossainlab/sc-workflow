# 04_differential_expression.R
# ==============================
# Pseudobulk differential expression and pathway enrichment.
#
# WHAT THIS SCRIPT DOES:
#   1. Cell proportion analysis — do cell type abundances shift between conditions?
#   2. Pseudobulk aggregation — sum raw counts per (sample x cell type)
#   3. DESeq2 — negative binomial DE test per cell type
#   4. Volcano plot (PDF)
#   5. Pathway enrichment via enrichR (requires internet)
#
# WHY PSEUDOBULK?
#   Direct cell-level tests (Wilcoxon) treat each cell as independent,
#   inflating significance through pseudoreplication. Pseudobulk aggregates
#   cells by biological replicate (sample) first — statistically correct.
#   You need >= 3 samples per condition.
#
# HOW TO USE:
#   1. Edit CONFIGURATION below.
#   2. Run:  Rscript 04_differential_expression.R
#
# OUTPUT:
#   results/04_de/de_results.csv
#   results/04_de/pathway_enrichment.csv
#   results/04_de/plots/

suppressMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(enrichR)
})


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH <- "results/03_annotation/seurat_annotated.rds"

# Column names in Seurat metadata
CELL_TYPE_COL <- "cell_type"
SAMPLE_COL    <- "sample"
CONDITION_COL <- "condition"

# --- Comparison ---------------------------------------------------------------
REFERENCE_CONDITION <- "Control"    # denominator (control / untreated)
TEST_CONDITION      <- "Treatment"  # numerator  (treated)

# Leave empty ("") to run DE for all cell types.
TARGET_CELL_TYPE <- ""

# --- DE thresholds ------------------------------------------------------------
PADJ_CUTOFF          <- 0.05
LFC_CUTOFF           <- 0.5
MIN_PSEUDOBULK_CELLS <- 10L   # minimum cells per sample to include in pseudobulk

# --- Pathway enrichment -------------------------------------------------------
RUN_ENRICHR  <- TRUE
ENRICHR_LIBS <- c("MSigDB_Hallmark_2020", "KEGG_2021_Human")

# --- Output -------------------------------------------------------------------
OUTPUT_DIR <- "results/04_de"


# ==============================================================================
#  STEP 1 — LOAD & VALIDATE
# ==============================================================================

load_data <- function() {
  if (!file.exists(INPUT_PATH)) stop("File not found: ", INPUT_PATH)
  seu <- readRDS(INPUT_PATH)

  required <- c(CELL_TYPE_COL, SAMPLE_COL, CONDITION_COL)
  missing  <- setdiff(required, colnames(seu@meta.data))
  if (length(missing) > 0L) {
    stop(sprintf(
      "Missing metadata columns: %s\nAvailable: %s",
      paste(missing, collapse = ", "),
      paste(colnames(seu@meta.data), collapse = ", ")
    ))
  }
  seu
}


# ==============================================================================
#  STEP 2 — CELL PROPORTION ANALYSIS
# ==============================================================================

cell_proportion_analysis <- function(seu, plot_dir) {
  meta <- seu@meta.data[, c(CELL_TYPE_COL, SAMPLE_COL, CONDITION_COL)]
  colnames(meta) <- c("cell_type", "sample", "condition")

  props <- meta |>
    count(sample, cell_type, condition) |>
    group_by(sample) |>
    mutate(proportion = n / sum(n)) |>
    ungroup()

  p <- ggplot(props, aes(sample, proportion, fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~condition, scales = "free_x") +
    labs(title = "Cell Type Proportions per Sample",
         y = "Fraction of cells", x = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(plot_dir, "01_cell_proportions.pdf"), p, width = 10, height = 5)
}


# ==============================================================================
#  STEP 3 — PSEUDOBULK DE (DESeq2)
# ==============================================================================

pseudobulk_de <- function(seu) {
  all_types    <- unique(as.character(seu@meta.data[[CELL_TYPE_COL]]))
  types_to_run <- if (nchar(TARGET_CELL_TYPE) > 0L) TARGET_CELL_TYPE else all_types

  all_results <- list()

  for (ct in types_to_run) {
    ct_mask <- seu@meta.data[[CELL_TYPE_COL]] == ct
    sub     <- seu[, ct_mask]
    meta    <- sub@meta.data

    # Pseudobulk: sum raw counts per sample
    counts_mat <- GetAssayData(sub, layer = "counts")
    samples     <- unique(as.character(meta[[SAMPLE_COL]]))

    pb_counts <- list()
    pb_meta   <- list()

    for (s in samples) {
      s_mask <- as.character(meta[[SAMPLE_COL]]) == s
      if (sum(s_mask) < MIN_PSEUDOBULK_CELLS) next
      pb_counts[[s]] <- Matrix::rowSums(counts_mat[, s_mask, drop = FALSE])
      pb_meta[[s]]   <- data.frame(
        sample    = s,
        condition = unique(as.character(meta[[CONDITION_COL]][s_mask]))[1]
      )
    }

    if (length(pb_counts) < 2L) next

    count_df <- do.call(cbind, pb_counts)
    meta_df  <- do.call(rbind, pb_meta)
    rownames(meta_df) <- meta_df$sample

    # Only keep samples from the two conditions being compared
    keep_samples <- meta_df$condition %in% c(REFERENCE_CONDITION, TEST_CONDITION)
    count_df <- count_df[, keep_samples, drop = FALSE]
    meta_df  <- meta_df[keep_samples, , drop = FALSE]

    if (length(unique(meta_df$condition)) < 2L) next

    # Remove genes with no counts
    count_df <- count_df[rowSums(count_df) >= 10L, , drop = FALSE]
    if (nrow(count_df) < 10L) next

    meta_df$condition <- relevel(factor(meta_df$condition), ref = REFERENCE_CONDITION)

    res <- tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = round(count_df),
        colData   = meta_df,
        design    = ~condition
      )
      dds <- DESeq(dds, quiet = TRUE)
      res <- results(dds,
                     contrast = c("condition", TEST_CONDITION, REFERENCE_CONDITION),
                     alpha    = PADJ_CUTOFF)
      df <- as.data.frame(res)
      df$gene      <- rownames(df)
      df$cell_type <- ct
      df
    }, error = function(e) {
      message(sprintf("  %s: DE failed — %s", ct, conditionMessage(e)))
      NULL
    })

    if (!is.null(res)) all_results[[ct]] <- res
  }

  if (length(all_results) == 0L) {
    stop("No DE results. Check REFERENCE_CONDITION/TEST_CONDITION match metadata values.")
  }

  de_results <- bind_rows(all_results)
  de_results$significant <- !is.na(de_results$padj) &
    de_results$padj <= PADJ_CUTOFF &
    abs(de_results$log2FoldChange) >= LFC_CUTOFF
  de_results
}


# ==============================================================================
#  STEP 4 — VOLCANO PLOT
# ==============================================================================

save_volcano_plot <- function(de_results, plot_dir) {
  df <- de_results |>
    filter(!is.na(padj)) |>
    mutate(
      neg_log10_padj = pmin(-log10(padj), 300),
      colour = case_when(
        log2FoldChange >  LFC_CUTOFF & padj <= PADJ_CUTOFF ~ "Up",
        log2FoldChange < -LFC_CUTOFF & padj <= PADJ_CUTOFF ~ "Down",
        TRUE ~ "NS"
      )
    )

  top_genes <- df |>
    filter(significant) |>
    group_by(cell_type) |>
    slice_min(padj, n = 10L) |>
    ungroup()

  p <- ggplot(df, aes(log2FoldChange, neg_log10_padj, colour = colour)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_text_repel(data = top_genes, aes(label = gene), size = 2.5,
                    max.overlaps = 15, show.legend = FALSE) +
    scale_colour_manual(values = c("Up" = "#C44E52", "Down" = "#4C72B0", "NS" = "#AAAAAA")) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", colour = "grey50") +
    labs(title = sprintf("Volcano: %s vs %s", TEST_CONDITION, REFERENCE_CONDITION),
         x = "log2 Fold Change", y = "-log10(padj)", colour = NULL) +
    facet_wrap(~cell_type) +
    theme_classic()

  n_types <- length(unique(df$cell_type))
  ncols   <- ceiling(sqrt(n_types))
  ggsave(file.path(plot_dir, "02_volcano.pdf"), p,
         width = 5 * ncols, height = 4 * ceiling(n_types / ncols))
}


# ==============================================================================
#  STEP 5 — PATHWAY ENRICHMENT
# ==============================================================================

run_pathway_enrichment <- function(de_results) {
  if (!RUN_ENRICHR) return(NULL)

  up_genes <- de_results$gene[!is.na(de_results$padj) &
    de_results$log2FoldChange >  LFC_CUTOFF &
    de_results$padj           <= PADJ_CUTOFF]
  dn_genes <- de_results$gene[!is.na(de_results$padj) &
    de_results$log2FoldChange < -LFC_CUTOFF &
    de_results$padj           <= PADJ_CUTOFF]

  collect <- list()
  for (direction in c("Up", "Down")) {
    genes <- if (direction == "Up") up_genes else dn_genes
    if (length(genes) < 5L) next
    tryCatch({
      enr <- enrichr(genes, ENRICHR_LIBS)
      for (lib in names(enr)) {
        r             <- head(enr[[lib]], 10L)
        r$direction   <- direction
        r$library     <- lib
        collect[[paste(direction, lib)]] <- r
      }
    }, error = function(e) {
      message(sprintf("  Enrichr (%s) failed: %s", direction, conditionMessage(e)))
    })
  }
  if (length(collect) == 0L) return(NULL)
  bind_rows(collect)
}


# ==============================================================================
#  MAIN
# ==============================================================================

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(OUTPUT_DIR, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  seu <- load_data()
  cell_proportion_analysis(seu, plot_dir)
  de_results <- pseudobulk_de(seu)
  save_volcano_plot(de_results, plot_dir)

  # Print top significant genes
  top_sig <- de_results |>
    filter(significant) |>
    arrange(padj) |>
    select(gene, cell_type, log2FoldChange, padj) |>
    head(20L)
  if (nrow(top_sig) > 0L) print(top_sig)

  enrichment_df <- run_pathway_enrichment(de_results)

  write.csv(de_results,    file.path(OUTPUT_DIR, "de_results.csv"),         row.names = FALSE)
  if (!is.null(enrichment_df)) {
    write.csv(enrichment_df, file.path(OUTPUT_DIR, "pathway_enrichment.csv"), row.names = FALSE)
  }
  message("Saved -> ", OUTPUT_DIR, "/")
}

main()
