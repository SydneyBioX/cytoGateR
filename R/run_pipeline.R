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
#' Build a marker-based, hierarchical soft-gating model for each lineage in a
#' marker table, compute per-cell membership probabilities, and derive hard
#' lineage/state labels.
#'
#' For each \code{cell_type} in \code{lineage_table}, a gating tree is fit from
#' positive and (optionally) negative marker sets. Per-cell probabilities are
#' then computed by evaluating each fitted tree via \code{tree_prob()}, with an
#' optional negative-marker penalty controlled by \code{neg_strength}.
#'
#' Hard labels are assigned as the lineage with the highest probability (excluding
#' \code{"Proliferating"}), unless the best probability is below
#' \code{uncert_thresh}, in which case the cell is labeled \code{"Uncertain"}.
#' If a \code{"Proliferating"} lineage is present, cells with proliferating
#' probability \eqn{\ge 0.5} and not \code{"Uncertain"} receive a \code{"Prolif "}
#' prefix in \code{cell_type_state}.
#'
#' @param spe A \code{SpatialExperiment} or \code{SingleCellExperiment} containing
#'   per-cell marker expression.
#' @param lineage_table A \code{data.frame}/tibble with at least the columns
#'   \code{cell_type}, \code{pos_markers}, and \code{neg_markers}. Marker columns
#'   must be list-columns of character vectors. Markers not present in
#'   \code{spe} are dropped during cleaning.
#' @param assay_name Character scalar naming the assay in \code{spe} used as the
#'   expression matrix (default \code{"norm"}).
#' @param max_depth Integer; maximum recursion depth of each gating tree.
#' @param min_cells Integer; minimum number of cells required to allow a node split.
#' @param min_score Numeric; minimum separability score required to accept a split.
#' @param uncert_thresh Numeric in [0,1]; cells whose best (non-proliferating)
#'   lineage probability is below this threshold are labeled \code{"Uncertain"}.
#' @param neg_strength Numeric in [0,1]; strength of the negative-marker penalty
#'   applied within \code{tree_prob()} (default \code{0.4}). Larger values increase
#'   the penalty from negative-marker expression.
#' @param cutoff_method Cutoff selection method passed to \code{fit_gmm_2()}.
#'   One of \code{"mean"} or \code{"equal_posteriors"}.
#' @param gmm_model_names Optional character vector of model names to pass to
#'   \code{mclust::Mclust()} for 1D GMM fitting (e.g., \code{"V"}, \code{"E"}).
#'   Use \code{NULL} (default) to allow selection by BIC.
#' @param workers Optional integer number of workers for parallel execution.
#'   Defaults to 1.
#'
#' @return A named list with components:
#' \describe{
#'   \item{spe}{The input \code{spe} with added \code{colData} columns:
#'     \code{P_*} (per-lineage probabilities), \code{cell_type_hard},
#'     \code{cell_type_state}, \code{prob_best}, and \code{prob_prolif}.}
#'   \item{trees}{A named list of fitted gating trees, one per \code{cell_type}.}
#'   \item{prob_mat}{Numeric matrix of per-cell probabilities
#'     (cells \eqn{\times} cell types).}
#'   \item{hard_label}{Character vector of hard lineage labels after applying
#'     \code{uncert_thresh}.}
#'   \item{lineage_table}{The cleaned marker table actually used to build the trees.}
#' }
#'
#' @details
#' Marker statistics used for logistic scoring/penalization are precomputed once
#' across all markers appearing in \code{lineage_table} via \code{fit_marker_stats()}.
#' Trees are built with \code{build_fullcoverage_tree()} and evaluated with
#' \code{tree_prob()}.
#'
#' @export
run_tree_gating <- function(spe,
                            lineage_table,
                            assay_name = "norm",
                            max_depth = 4,
                            min_cells = 200,
                            min_score = 0.5,
                            uncert_thresh = 0.25,
                            neg_strength = 0.4, # <--- ADD HERE
                            cutoff_method = c("mean", "equal_posteriors"),
                            gmm_model_names = NULL,
                            # parallel = FALSE,
                            workers = 1) {
  cutoff_method <- match.arg(cutoff_method)

  .assert_spe(spe)
  .assert_lineage_table(lineage_table)

  lineage_table <- .clean_lineage_table(lineage_table, spe)

  # --- ADD THIS LINE ---
  # Pre-calculate stats for all markers used in the lineage table
  all_markers <- unique(unlist(c(lineage_table$pos_markers, lineage_table$neg_markers))) #
  marker_stats <- fit_marker_stats(spe, all_markers, assay_name = assay_name) #


  expr_norm <- SummarizedExperiment::assay(spe, assay_name)

  build_tree <- function(i) {
    ct <- lineage_table$cell_type[i]
    pos <- lineage_table$pos_markers[[i]]
    neg <- lineage_table$neg_markers[[i]] %||% character(0)
    message(paste0("Building ", i, "th tree - ", ct, " Cell Type."))
    timing <- system.time(
      tree <- build_fullcoverage_tree(
        expr_norm,
        pos,
        neg,
        max_depth = max_depth,
        min_cells = min_cells,
        min_score = min_score,
        cutoff_method = cutoff_method,
        gmm_model_names = gmm_model_names
      )
    )
    message("Finished ",
            ct,
            ". Timing: ",
            paste(
              c("usr", "sys", "elap", "u.c", "s.c"),
              round(timing, 2),
              sep = ": ",
              collapse = "; "
            ))
    tree
  }

  if (workers > 1) {
    oplan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(oplan), add = TRUE)
    trees <- furrr::future_map(
      seq_len(nrow(lineage_table)),
      build_tree,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    trees <- lapply(seq_len(nrow(lineage_table)), build_tree)
  }

  trees <- setNames(trees, lineage_table$cell_type)

  prob_mat <- sapply(names(trees), function(ct) {
    neg <- lineage_table$neg_markers[lineage_table$cell_type == ct][[1]] %||% character(0)
    vapply(seq_len(ncol(expr_norm)), function(i) {
      tree_prob(
        trees[[ct]],
        expr_norm,
        i,
        combine = "mean",
        neg_markers = neg,
        marker_stats = marker_stats,
        neg_strength = neg_strength # <--- PASS HERE
      )
    }, numeric(1))
  })

  prob_mat[!is.finite(prob_mat)] <- 0

  lineage_priority <- lineage_table$cell_type[lineage_table$cell_type != "Proliferating"]
  prob_lineage <- prob_mat[, lineage_priority, drop = FALSE]

  if (ncol(prob_lineage) == 0) {
    best_lab <- rep("Uncertain", nrow(prob_mat))
    best_p <- rep(0, nrow(prob_mat))
  } else {
    best_idx <- apply(prob_lineage, 1, which.max)
    best_lab <- colnames(prob_lineage)[best_idx]
    best_p <- apply(prob_lineage, 1, max)
  }

  hard_label <- ifelse(best_p < uncert_thresh, "Uncertain", best_lab)

  if ("Proliferating" %in% colnames(prob_mat)) {
    p_prolif <- prob_mat[, "Proliferating"]
  } else {
    p_prolif <- rep(0, nrow(prob_mat))
  }
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
