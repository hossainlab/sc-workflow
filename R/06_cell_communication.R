# 06_cell_communication.R
# =========================
# Infer intercellular signalling using CellChat.
#
# WHAT THIS SCRIPT DOES:
#   1. Create a CellChat object from the annotated Seurat object
#   2. Compute communication probabilities (all ligand-receptor pairs)
#   3. Infer signalling pathways
#   4. Aggregate network and identify dominant senders / receivers
#   5. (Optional) Run CellChat per condition for comparative analysis
#   6. Save results
#
# IMPORTANT CAVEAT:
#   CellChat predicts *potential* communication from co-expression of ligands
#   and receptors. Results are hypotheses that require wet-lab validation.
#
# HOW TO USE:
#   1. Edit CONFIGURATION below.
#   2. Run:  Rscript 06_cell_communication.R
#
# OUTPUT:
#   results/06_communication/cellchat.rds
#   results/06_communication/interactions.csv
#   results/06_communication/plots/

suppressMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)


# ==============================================================================
#  CONFIGURATION — Edit this section before running
# ==============================================================================

INPUT_PATH    <- "results/03_annotation/seurat_annotated.rds"
CELL_TYPE_COL <- "cell_type"

# Column for condition — leave "" for single-condition analysis.
CONDITION_COL <- "condition"

# --- CellChat parameters ------------------------------------------------------
# CellChat database. Options: "CellChatDB.human" or "CellChatDB.mouse"
CELLCHAT_DB <- "CellChatDB.human"

# Minimum fraction of cells in a group that must express the ligand/receptor.
MIN_EXPR_PROP <- 0.1

# Number of CPU cores for permutation testing.
N_WORKERS <- 4L

# Number of top interactions to display.
N_TOP_INTERACTIONS <- 20L

# --- Comparative analysis -----------------------------------------------------
# Set TRUE to run CellChat separately per condition (requires CONDITION_COL).
RUN_COMPARATIVE <- FALSE

# --- Output -------------------------------------------------------------------
OUTPUT_DIR <- "results/06_communication"


# ==============================================================================
#  STEP 1 — LOAD DATA
# ==============================================================================

load_data <- function() {
  if (!file.exists(INPUT_PATH)) stop("File not found: ", INPUT_PATH)
  seu <- readRDS(INPUT_PATH)
  if (!CELL_TYPE_COL %in% colnames(seu@meta.data)) {
    stop(sprintf("Column '%s' not found in metadata.", CELL_TYPE_COL))
  }
  seu
}


# ==============================================================================
#  STEP 2 — CREATE CELLCHAT OBJECT
# ==============================================================================

create_cellchat <- function(seu) {
  db    <- get(CELLCHAT_DB, envir = asNamespace("CellChat"))
  interaction_input <- GetAssayData(seu, layer = "data")
  meta  <- data.frame(group = as.character(seu@meta.data[[CELL_TYPE_COL]]),
                      row.names = colnames(seu))

  cc <- createCellChat(object = interaction_input, meta = meta, group.by = "group")
  CellChatDB(cc) <- subsetDB(db, search = "Secreted Signaling")
  cc
}


# ==============================================================================
#  STEP 3 — COMPUTE COMMUNICATION PROBABILITIES
# ==============================================================================

run_cellchat <- function(cc) {
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, thresh.pc = MIN_EXPR_PROP)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10L)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  cc
}


# ==============================================================================
#  STEP 4 — VISUALISATION
# ==============================================================================

save_cellchat_plots <- function(cc, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Interaction number heatmap
  pdf(file.path(plot_dir, "01_interaction_heatmap.pdf"), width = 8, height = 7)
  netVisual_heatmap(cc, measure = "count", color.heatmap = "Reds")
  dev.off()

  # Interaction weight/strength heatmap
  pdf(file.path(plot_dir, "02_weight_heatmap.pdf"), width = 8, height = 7)
  netVisual_heatmap(cc, measure = "weight", color.heatmap = "GnBu")
  dev.off()

  # Sender / receiver dominance scatter
  p_sc <- netAnalysis_signalingRole_scatter(cc)
  ggsave(file.path(plot_dir, "03_sender_receiver.pdf"), p_sc, width = 7, height = 6)

  # Top signalling pathways bubble chart
  p_bub <- netVisual_bubble(cc, remove.isolate = FALSE,
                             thresh = 0.05, n.top = N_TOP_INTERACTIONS)
  ggsave(file.path(plot_dir, "04_top_interactions_bubble.pdf"), p_bub,
         width = 10, height = 6)
}


# ==============================================================================
#  STEP 5 — COMPARATIVE ANALYSIS (OPTIONAL)
# ==============================================================================

run_comparative_cellchat <- function(seu) {
  if (!RUN_COMPARATIVE) return(NULL)
  if (nchar(CONDITION_COL) == 0L || !CONDITION_COL %in% colnames(seu@meta.data)) {
    message("Condition column not found — skipping comparative analysis.")
    return(NULL)
  }

  conditions <- unique(as.character(seu@meta.data[[CONDITION_COL]]))
  cc_list    <- list()

  for (cond in conditions) {
    mask <- as.character(seu@meta.data[[CONDITION_COL]]) == cond
    sub  <- seu[, mask]
    cc   <- create_cellchat(sub)
    cc   <- run_cellchat(cc)
    cc_list[[cond]] <- cc
    message(sprintf("CellChat complete for condition: %s", cond))
  }

  combined <- mergeCellChat(cc_list, add.names = conditions)
  combined
}


# ==============================================================================
#  MAIN
# ==============================================================================

main <- function() {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(OUTPUT_DIR, "plots")

  seu <- load_data()
  cc  <- create_cellchat(seu)
  cc  <- run_cellchat(cc)
  save_cellchat_plots(cc, plot_dir)

  # Save flat interaction table
  interactions <- subsetCommunication(cc)
  write.csv(interactions, file.path(OUTPUT_DIR, "interactions.csv"), row.names = FALSE)

  # Print top interactions
  top_ix <- head(interactions[order(interactions$prob, decreasing = TRUE), ], N_TOP_INTERACTIONS)
  message("\nTop interactions:")
  print(top_ix[, c("source", "target", "ligand", "receptor", "prob", "pathway_name")])

  saveRDS(cc, file.path(OUTPUT_DIR, "cellchat.rds"))

  # Comparative analysis
  combined <- run_comparative_cellchat(seu)
  if (!is.null(combined)) {
    saveRDS(combined, file.path(OUTPUT_DIR, "cellchat_combined.rds"))
  }

  message("Saved -> ", OUTPUT_DIR, "/")
}

main()
