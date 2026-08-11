#!/usr/bin/env Rscript

# Export R->Python parity fixtures from the stored tree-gating reference objects.
#
# The reference RDS files under
#   PRJ-cytogater/20260624-python_pacakge_tests-LY/in/tree_gated_reference/
# are the saved return value of cytoGateR::run_tree_gating(), produced by
#   20260714-constract_reference-LY/src/06-construct_reference-ML.R
# with the parameters recorded in R_PARAMS below.
#
# This script only reads those objects, so it needs SingleCellExperiment but NOT
# cytoGateR itself -- which matters because cytoGateR is not installed on every
# node. Re-run it whenever the reference objects are regenerated.
#
# Usage:
#   Rscript tests/fixtures/export_r_reference.R [dataset ...]

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(jsonlite)
})

REFERENCE_DIR <- "/dskh/nobackup/biostat/projects/PRJ-cytogater/20260624-python_pacakge_tests-LY/in/tree_gated_reference"
EXPRESSION_DIR <- "/dskh/nobackup/biostat/projects/PRJ-cytogater/20260624-python_pacakge_tests-LY/in/datasets_csv"
OUT_DIR <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "r_parity")

# The exact run_tree_gating() call the reference objects were built with.
R_PARAMS <- list(
  assay_name = "chosen_assay",
  max_depth = 6L,
  min_cells = 100L,
  min_score = 0.05,
  uncert_thresh = 0.25,
  neg_strength = 0.5,
  cutoff_method = "equal_posteriors",
  gmm_model_names = NULL
)

# Number of cells retained in the row-level fixtures. The Python side still runs
# on the full matrix -- the GMM cutoffs are fit across all cells, so subsetting
# the *input* would change the answer. Only the comparison is subsetted.
N_SAMPLE <- 500L

datasets <- commandArgs(trailingOnly = TRUE)
if (length(datasets) == 0) datasets <- c("breast_MIBI-TOF_sce")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (dataset in datasets) {
  message(sprintf("=== %s ===", dataset))
  reference_path <- file.path(REFERENCE_DIR, paste0(dataset, ".rds"))
  expression_path <- file.path(EXPRESSION_DIR, paste0(dataset, ".csv"))
  stopifnot(file.exists(reference_path), file.exists(expression_path))

  result <- readRDS(reference_path)
  spe <- result$spe
  prob_mat <- result$prob_mat
  hard_label <- result$hard_label
  cell_ids <- colnames(spe)
  cell_types <- colnames(prob_mat)

  # The Python test reads expression from the shared CSV, so confirm the CSV is
  # byte-identical to the assay the reference was gated on before trusting it.
  assay_matrix <- assay(spe, R_PARAMS$assay_name)
  csv_head <- read.csv(expression_path, nrows = 200, check.names = FALSE)
  csv_markers <- setdiff(names(csv_head), "cell_id")
  shared_markers <- intersect(csv_markers, rownames(assay_matrix))
  csv_values <- as.matrix(csv_head[, shared_markers, drop = FALSE])
  assay_values <- t(assay_matrix[shared_markers, as.character(csv_head$cell_id), drop = FALSE])
  max_diff <- max(abs(csv_values - assay_values))
  message(sprintf("  CSV vs assay max abs diff (first 200 cells): %.3g", max_diff))
  stopifnot(max_diff < 1e-9)

  # Evenly spaced rather than random so the fixture is reproducible without a seed.
  keep <- unique(round(seq(1, length(cell_ids), length.out = min(N_SAMPLE, length(cell_ids)))))

  probs_out <- data.frame(cell_id = cell_ids[keep], stringsAsFactors = FALSE)
  probs_out <- cbind(probs_out, as.data.frame(prob_mat[keep, , drop = FALSE]))
  write.csv(probs_out, file.path(OUT_DIR, paste0(dataset, "__probs_sample.csv")), row.names = FALSE)

  write.csv(
    data.frame(cell_id = cell_ids[keep], hard_label = hard_label[keep], stringsAsFactors = FALSE),
    file.path(OUT_DIR, paste0(dataset, "__labels_sample.csv")),
    row.names = FALSE
  )

  lineage <- result$lineage_table
  # I() keeps length-1 marker vectors as JSON arrays. Without it auto_unbox
  # turns them into bare strings, and the Python side then iterates the string
  # character by character and silently drops the cell type.
  lineage_out <- lapply(seq_len(nrow(lineage)), function(i) {
    list(
      cell_type = as.character(lineage$cell_type[[i]]),
      pos_markers = I(as.character(lineage$pos_markers[[i]])),
      neg_markers = I(as.character(lineage$neg_markers[[i]]))
    )
  })
  writeLines(
    toJSON(lineage_out, auto_unbox = TRUE, pretty = TRUE),
    file.path(OUT_DIR, paste0(dataset, "__lineage.json"))
  )

  # Whole-dataset aggregates catch drift that a 500-cell sample would miss.
  label_counts <- as.list(table(hard_label))
  summary <- list(
    dataset = dataset,
    expression_csv = expression_path,
    reference_rds = reference_path,
    params = R_PARAMS,
    n_cells = ncol(spe),
    n_markers = nrow(spe),
    cell_types = cell_types,
    sample_indices = keep,
    label_counts = label_counts,
    prob_col_means = as.list(setNames(colMeans(prob_mat), cell_types)),
    prob_col_sds = as.list(setNames(apply(prob_mat, 2, sd), cell_types))
  )
  writeLines(
    toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = 15, null = "null"),
    file.path(OUT_DIR, paste0(dataset, "__summary.json"))
  )

  message(sprintf("  wrote fixtures for %d cells x %d cell types (%d sampled)",
                  ncol(spe), length(cell_types), length(keep)))
}

message("done")
