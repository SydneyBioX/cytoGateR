#' Run soft-gating pipeline end-to-end
#'
#' @param spe SpatialExperiment / SingleCellExperiment.
#' @param lineage_table Tibble with cell_type/pos_markers/neg_markers.
#' @param assay_name Assay to use (default "norm").
#' @param unknown_thresh Threshold for Unknown label (default 0.4).
#' @param store Store results back into spe (default TRUE).
#'
#' @return A list with: spe (if store=TRUE), marker_stats, prob_mat, labels.
#' @export
run_soft_gating <- function(spe, lineage_table,
                            assay_name = "norm",
                            unknown_thresh = 0.4,
                            store = TRUE) {
  .assert_spe(spe)
  .assert_lineage_table(lineage_table)

  # filter markers to those present
  # lineage_table2 <- lineage_table |>
  #   dplyr::mutate(
  #     pos_markers = lapply(pos_markers, function(v) intersect(unlist(v), rownames(spe))),
  #     neg_markers = lapply(neg_markers, function(v) intersect(unlist(v), rownames(spe)))
  #   ) |>
  #   dplyr::filter(lengths(pos_markers) > 0)

  lineage_table2 <- .clean_lineage_table(lineage_table, spe)




  all_markers <- unique(unlist(c(lineage_table2$pos_markers, lineage_table2$neg_markers)))

  marker_stats <- fit_marker_stats(spe, all_markers, assay_name = assay_name)
  prob_mat <- calculate_soft_scores(spe, marker_stats, lineage_table2, assay_name = assay_name)
  prob_mat[!is.finite(prob_mat)] <- 0

  labels <- assign_soft_labels(prob_mat, unknown_thresh = unknown_thresh)

  if (store) {
    # Store prob_mat as a SingleCellExperiment assay (recommended)
    for (ct in colnames(prob_mat)) {
      SummarizedExperiment::colData(spe)[[paste0("P_", gsub("\\s+", "_", ct))]] <- prob_mat[, ct]
    }
    SummarizedExperiment::colData(spe)$soft_tree_label <- labels


    # Store labels in colData
    SummarizedExperiment::colData(spe)$soft_tree_label <- labels
  }

  list(
    spe = spe,
    marker_stats = marker_stats,
    prob_mat = prob_mat,
    labels = labels,
    lineage_table = lineage_table2
  )
}


#' Run tree-based soft gating on a SpatialExperiment
#'
#' Fit marker-based gating trees for multiple cell types and compute per-cell
#' membership probabilities and hard labels using a hierarchical, soft
#' decision-tree model.
#'
#' @param spe A `SpatialExperiment` or `SingleCellExperiment` containing per-cell
#'   marker expression.
#' @param lineage_table A tibble/data.frame with columns `cell_type`,
#'   `pos_markers`, and `neg_markers`. Marker columns must be list-columns of
#'   character vectors specifying positive and negative markers for each cell
#'   type.
#' @param assay_name Name of the assay in `spe` to use as the expression matrix
#'   (e.g. `"norm"`).
#' @param max_depth Maximum recursion depth of each gating tree.
#' @param min_cells Minimum number of cells required to allow a node split.
#' @param min_score Minimum separability score required to accept a split.
#' @param uncert_thresh Cells whose best lineage probability is below this
#'   threshold are labeled `"Uncertain"`.
#'
#' @return A named list with components:
#' \describe{
#'   \item{spe}{The input `spe` object with additional columns added to `colData`
#'     (`P_*` probabilities, `cell_type_hard`, `cell_type_state`, `prob_best`,
#'     `prob_prolif`).}
#'   \item{trees}{A named list of fitted gating trees, one per cell type.}
#'   \item{prob_mat}{Numeric matrix of per-cell probabilities
#'     (cells × cell types).}
#'   \item{hard_label}{Character vector of hard lineage labels after applying the
#'     uncertainty threshold.}
#'   \item{lineage_table}{The cleaned marker table actually used to build the
#'     trees (after intersecting with available markers in `spe`).}
#' }
#'
#' @export
run_tree_gating <- function(spe,
                            lineage_table,
                            assay_name = "norm",
                            max_depth = 4,
                            min_cells = 200,
                            min_score = 0.5,
                            uncert_thresh = 0.25) {

  .assert_spe(spe)
  .assert_lineage_table(lineage_table)

  lineage_table <- .clean_lineage_table(lineage_table, spe)

  expr_norm <- SummarizedExperiment::assay(spe, assay_name)

  # Setup parallel processing with future
  oplan <- future::plan(future::multisession, workers = parallel::detectCores() - 1)
  on.exit(future::plan(oplan), add = TRUE)

  # Build trees in parallel
  trees_list <- furrr::future_map(
    seq_len(nrow(lineage_table)),
    function(i) {
      pos <- lineage_table$pos_markers[[i]]
      neg <- lineage_table$neg_markers[[i]] %||% character(0)

      build_fullcoverage_tree(expr_norm, pos, neg,
                              max_depth = max_depth,
                              min_cells = min_cells,
                              min_score = min_score)
    },
    .options = furrr::furrr_options(seed = TRUE)
  )

  trees <- setNames(trees_list, lineage_table$cell_type)

  prob_mat <- sapply(names(trees), function(ct) {
    neg <- lineage_table$neg_markers[lineage_table$cell_type == ct][[1]] %||% character(0)
    vapply(seq_len(ncol(expr_norm)), function(i) {
      tree_prob(trees[[ct]], expr_norm, i, combine = "mean", neg_markers = neg)
    }, numeric(1))
  })

  prob_mat[!is.finite(prob_mat)] <- 0

  lineage_priority <- setdiff(lineage_table$cell_type, "Proliferating")
  prob_lineage <- prob_mat[, lineage_priority, drop = FALSE]

  best_idx <- apply(prob_lineage, 1, which.max)
  best_lab <- colnames(prob_lineage)[best_idx]
  best_p <- apply(prob_lineage, 1, max)

  hard_label <- ifelse(best_p < uncert_thresh, "Uncertain", best_lab)

  has_prolif <- "Proliferating" %in% colnames(prob_mat)
  p_prolif <- if (has_prolif) prob_mat[, "Proliferating"] else rep(0, nrow(prob_mat))
  prolif_flag <- p_prolif >= 0.5

  hard_label_with_state <- ifelse(prolif_flag & hard_label != "Uncertain",
                                  paste("Prolif", hard_label),
                                  hard_label)

  SummarizedExperiment::colData(spe)$prob_best <- best_p
  SummarizedExperiment::colData(spe)$cell_type_hard <- factor(hard_label)
  SummarizedExperiment::colData(spe)$cell_type_state <- factor(hard_label_with_state)
  SummarizedExperiment::colData(spe)$prob_prolif <- p_prolif

  for (ct in colnames(prob_mat)) {
    SummarizedExperiment::colData(spe)[[paste0("P_", gsub("\\s+", "_", ct))]] <- prob_mat[, ct]
  }

  list(spe = spe, trees = trees, prob_mat = prob_mat, hard_label = hard_label,
       lineage_table = lineage_table)
}
