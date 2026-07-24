#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @keywords internal
.safe_mad <- function(x, na.rm = TRUE) {
  v <- stats::mad(x, na.rm = na.rm)
  if (!is.finite(v) || v == 0) 1e-6 else v
}

#' @keywords internal
.assert_spe <- function(spe) {
  ok <- inherits(spe, c("SpatialExperiment", "SingleCellExperiment"))
  if (!ok) stop("spe must be a SpatialExperiment or SingleCellExperiment.")
}

#' @keywords internal
.assert_lineage_table <- function(lineage_table) {
  req <- c("cell_type", "pos_markers", "neg_markers")
  miss <- setdiff(req, names(lineage_table))
  if (length(miss) > 0) stop("lineage_table missing columns: ", paste(miss, collapse = ", "))

  # minimal type checks
  if (!is.list(lineage_table$pos_markers) || !is.list(lineage_table$neg_markers)) {
    stop("lineage_table$pos_markers and neg_markers must be list-columns (each element a character vector).")
  }
}

#' @keywords internal
.clean_lineage_table <- function(lineage_table, spe) {

  lineage_table <- dplyr::as_tibble(lineage_table)

  lineage_table |>
    dplyr::mutate(
      cell_type   = trimws(as.character(.data$cell_type)),
      pos_markers = lapply(.data$pos_markers, function(v) intersect(unlist(v), rownames(spe))),
      neg_markers = lapply(.data$neg_markers, function(v) intersect(unlist(v), rownames(spe)))
    ) |>
    dplyr::filter(lengths(.data$pos_markers) > 0)
}




#' Calculate F1 Score per Cell Type
#'
#' @param spe SpatialExperiment object
#' @param ref_col Column name for ground truth (e.g., "manual_label")
#' @param pred_col Column name for results (e.g., "hier_label")
#'
#' @return A data frame with Precision, Recall, and F1 score per type
#' @export
calculate_f1 <- function(spe, ref_col = "ref_broad", pred_col = "pred_broad") {
  df <- as.data.frame(SummarizedExperiment::colData(spe))

  # Get only the broad categories we care about
  types <- unique(as.character(df[[ref_col]]))
  types <- types[!is.na(types) & !grepl("unassigned|undefined|unknown", types, ignore.case = TRUE)]

  results <- lapply(types, function(type) {
    tp <- sum(df[[pred_col]] == type & df[[ref_col]] == type, na.rm = TRUE)
    fp <- sum(df[[pred_col]] == type & df[[ref_col]] != type, na.rm = TRUE)
    fn <- sum(df[[pred_col]] != type & df[[ref_col]] == type, na.rm = TRUE)

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall    <- if ((tp + fn) > 0) tp / (tp + fn) else 0

    # It should be 2 * (p * r) / (p + r).
    f1_actual <- if ((precision + recall) > 0) (2 * precision * recall) / (precision + recall) else 0

    data.frame(
      Category = type,
      Precision = round(precision, 3),
      Recall = round(recall, 3),
      F1_Score = round(f1_actual, 3),
      Cell_Count = sum(df[[ref_col]] == type, na.rm = TRUE)
    )
  })

  do.call(rbind, results)
}






#' Calculate Full Uncertainty Suite with Sample Grouping
#'
#' @param prob_mat The [Cells x Types] matrix from your kNN or randomforest
#' @param spe The SpatialExperiment object.
#' @param sample_col The column identifying different images/samples.
#' @param k_spatial Neighbors for the spatial discordance check (default 15).
#'   from training) to rbind with prob_mat before processing. Both must share
#'   the same column names.
#' @export
calculate_uncertainty <- function(prob_mat,
                                  spe,
                                  sample_col = "sample_id",
                                  k_spatial = 15) {



  prob_mat <- as.matrix(prob_mat)

  # Reorder rows to match SPE cell order
  spe_cells <- colnames(spe)

  missing <- setdiff(spe_cells, rownames(prob_mat))
  extra   <- setdiff(rownames(prob_mat), spe_cells)

  if (length(missing) > 0) {
    warning(sprintf("%d cells in SPE are missing from prob_mat and will get NA uncertainty scores.", length(missing)))
  }
  if (length(extra) > 0) {
    warning(sprintf("%d rows in prob_mat have no matching SPE cell and will be dropped.", length(extra)))
  }

  # Align to SPE order — cells missing from prob_mat become NA rows
  common_cells  <- intersect(spe_cells, rownames(prob_mat))
  prob_aligned  <- matrix(NA,
                          nrow = length(spe_cells),
                          ncol = ncol(prob_mat),
                          dimnames = list(spe_cells, colnames(prob_mat)))
  prob_aligned[common_cells, ] <- prob_mat[common_cells, ]
  prob_mat <- prob_aligned

  # Entropy
  n_types <- ncol(prob_mat)
  labels  <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]

  # entropy <- -rowSums(prob_mat * log(prob_mat + 1e-10), na.rm = TRUE) / log(n_types)

  log_prob <- matrix(0, nrow = nrow(prob_mat), ncol = ncol(prob_mat))
  pos_mask <- prob_mat > 0
  log_prob[pos_mask] <- log(prob_mat[pos_mask])
  entropy <- -rowSums(prob_mat * log_prob, na.rm = TRUE) / log(n_types)


  gini_raw  <- 1 - rowSums(prob_mat^2, na.rm = TRUE)
  gini_norm <- gini_raw / (1 - 1/n_types)

  sorted_probs <- t(apply(prob_mat, 1, sort, decreasing = TRUE))
  margin_val   <- 1 - (sorted_probs[, 1] - sorted_probs[, 2])

  # Spatial Metrics (Sample-Aware)
  samples        <- SummarizedExperiment::colData(spe)[[sample_col]]
  unique_samples <- unique(samples)
  spatial_discordance <- rep(NA, nrow(prob_mat))

  for (s in unique_samples) {
    idx <- which(samples == s)
    if (length(idx) <= k_spatial) next

    coords_subset <- SpatialExperiment::spatialCoords(spe)[idx, ]
    knn_res       <- dbscan::kNN(coords_subset, k = k_spatial)
    sample_labels <- labels[idx]

    sample_discordance <- sapply(seq_len(length(idx)), function(i) {
      neighbor_labels <- sample_labels[knn_res$id[i, ]]
      sum(neighbor_labels != sample_labels[i]) / k_spatial
    })
    spatial_discordance[idx] <- sample_discordance
  }


  data.frame(
    cell_id              = spe_cells,
    entropy              = entropy,
    gini_impurity        = gini_norm,
    margin_uncertainty   = margin_val,
    spatial_discordance  = spatial_discordance
  )
}





#' Apply Protected Spatial Priors to kNN Probabilities
#'
#' @param spe SpatialExperiment object.
#' @param prob_mat The full probability matrix [Cells x Types] from predict_unknown_with_knn.
#' @param k_spatial Number of physical neighbors for the "Prior" (default 50).
#' @param lambda The "Spatial Weight" (default 0.2).
#' @param protect_threshold Confidence level above which spatial priors are ignored (default 0.85).
#' @param out_col Column name for the new labels.
#' @export
calculate_spatial_prior_labels <- function(spe,
                                           prob_mat,
                                           k_spatial = 50,
                                           lambda = 0.2,
                                           protect_threshold = 0.85,
                                           out_col = "knn_spatial_label") {


  prob_mat  <- as.matrix(prob_mat)
  spe_cells <- colnames(spe)

  missing <- setdiff(spe_cells, rownames(prob_mat))
  extra   <- setdiff(rownames(prob_mat), spe_cells)

  if (length(missing) > 0) {
    warning(sprintf("%d cells in SPE are missing from prob_mat and will get NA labels.", length(missing)))
  }
  if (length(extra) > 0) {
    warning(sprintf("%d rows in prob_mat have no matching SPE cell and will be dropped.", length(extra)))
  }

  common_cells <- intersect(spe_cells, rownames(prob_mat))
  prob_aligned <- matrix(NA,
                         nrow = length(spe_cells),
                         ncol = ncol(prob_mat),
                         dimnames = list(spe_cells, colnames(prob_mat)))
  prob_aligned[common_cells, ] <- prob_mat[common_cells, ]
  prob_mat <- prob_aligned


  coords     <- SpatialExperiment::spatialCoords(spe)
  cell_types <- colnames(prob_mat)

  # Extract kNN Confidence (the max probability per row)
  knn_confidence      <- apply(prob_mat, 1, max)
  current_best_labels <- cell_types[max.col(prob_mat, ties.method = "first")]

  # Fast Spatial Neighbor Search
  knn_spatial <- dbscan::kNN(coords, k = k_spatial)

  message(sprintf("Integrating spatial context (Protecting cells > %s confidence)...", protect_threshold))

  # Apply Bayesian Update
  spatial_results <- vapply(seq_len(nrow(prob_mat)), function(i) {
    if (is.na(knn_confidence[i])) return(NA_character_)

    if (knn_confidence[i] >= protect_threshold) {
      return(current_best_labels[i])
    }

    neighbor_idx    <- knn_spatial$id[i, ]
    neighbor_labels <- current_best_labels[neighbor_idx]

    prior_counts <- table(factor(neighbor_labels, levels = cell_types))
    prior_probs  <- as.numeric(prior_counts) / k_spatial

    posterior <- prob_mat[i, ] * (prior_probs ^ lambda)
    posterior[!is.finite(posterior)] <- 0

    return(cell_types[which.max(posterior)])
  }, character(1))

  SummarizedExperiment::colData(spe)[[out_col]] <- spatial_results
  return(spe)
}
