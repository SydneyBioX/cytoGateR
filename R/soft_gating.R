#' Fit per-marker 2-GMM statistics for soft gating
#'
#' Fits a 2-component Gaussian Mixture Model (GMM) per marker and returns
#' a list of per-marker stats (cutoff, scale, weight, mu_high, valid).
#'
#' @param spe SpatialExperiment / SingleCellExperiment.
#' @param markers Character vector of marker names.
#' @param assay_name Assay to use (default "exprs").
#' @param min_n Minimum finite observations required to fit GMM (default 50).
#'
#' @return Named list. Each element is a list with fields:
#'   valid, cutoff, scale, weight, mu_high.
#'
#' @export
fit_marker_stats <- function(spe, markers, assay_name = "exprs", min_n = 50L) {
  .assert_spe(spe)
  if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
    stop("Assay '", assay_name, "' not found in spe.")
  }

  expr_mat <- SummarizedExperiment::assay(spe, assay_name)
  valid_markers <- intersect(markers, rownames(spe))

  stats_list <- vector("list", length(valid_markers))
  names(stats_list) <- valid_markers

  for (m in valid_markers) {
    vals <- as.numeric(expr_mat[m, ])
    vals <- vals[is.finite(vals)]

    # fallback defaults
    fallback <- list(
      valid  = FALSE,
      cutoff = stats::quantile(vals, 0.98, na.rm = TRUE),
      scale  = .safe_mad(vals, na.rm = TRUE),
      weight = 0.1,
      mu_high = if (length(vals) > 0) max(vals, na.rm = TRUE) else NA_real_
    )
    if (!is.finite(fallback$cutoff)) fallback$cutoff <- 0
    if (!is.finite(fallback$scale)  || fallback$scale == 0) fallback$scale <- 1e-3

    if (length(vals) < min_n || length(unique(vals)) < 3L) {
      stats_list[[m]] <- fallback
      next
    }

    gmm <- tryCatch(mclust::Mclust(vals, G = 2, verbose = FALSE), error = function(e) NULL)

    valid_fit <- FALSE
    if (!is.null(gmm) && length(gmm$parameters$mean) == 2) {
      mu <- as.numeric(gmm$parameters$mean)
      vars <- as.numeric(gmm$parameters$variance$sigmasq)

      # handle equal-variance models
      if (length(vars) == 1) vars <- rep(vars, 2)

      if (all(is.finite(mu)) && all(is.finite(vars)) && all(vars > 0)) {
        valid_fit <- TRUE
      }
    }

    if (!valid_fit) {
      stats_list[[m]] <- fallback
      next
    }

    ord <- order(mu)
    mu_low  <- mu[ord[1]]
    mu_high <- mu[ord[2]]
    sd_low  <- sqrt(vars[ord[1]])
    sd_high <- sqrt(vars[ord[2]])

    cutoff <- mean(c(mu_low, mu_high))
    scale_param <- mean(c(sd_low, sd_high))
    if (!is.finite(scale_param) || scale_param == 0) scale_param <- .safe_mad(vals)

    denom <- sqrt(sd_low^2 + sd_high^2)
    sep_score <- if (is.finite(denom) && denom > 1e-6) (mu_high - mu_low) / denom else 0
    if (!is.finite(sep_score) || sep_score <= 0) sep_score <- 0.1

    stats_list[[m]] <- list(
      valid  = TRUE,
      cutoff = cutoff,
      scale  = scale_param,
      weight = sep_score,
      mu_high = mu_high
    )
  }

  stats_list
}


#' Assign hard labels from probability matrix
#'
#' @param prob_mat Matrix [cells x types].
#' @param unknown_thresh If max prob < this, label as "Unknown" (default 0.4).
#'
#' @return Character vector of labels length nrow(prob_mat).
#' @export
assign_soft_labels <- function(prob_mat, unknown_thresh = 0.4) {
  if (!is.matrix(prob_mat)) stop("prob_mat must be a matrix.")
  if (ncol(prob_mat) == 0) stop("prob_mat has 0 columns.")

  max_probs <- apply(prob_mat, 1, max)
  lab <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
  lab[max_probs < unknown_thresh] <- "Unknown"
  lab
}



#' Calculate soft gating probabilities per cell type
#'
#' Computes a probability matrix [cells x cell_types] by combining
#' per-marker logistic scores weighted by separability, with optional
#' negative-marker penalties.
#'
#' @param spe SpatialExperiment / SingleCellExperiment.
#' @param marker_stats Output of fit_marker_stats().
#' @param lineage_table Tibble with columns: cell_type, pos_markers (list), neg_markers (list).
#' @param assay_name Assay to use (default "exprs").
#' @param neg_strength Penalty strength multiplier for negative markers (default 0.8).
#'
#' @return Numeric matrix [ncol(spe) x n_types].
#' @export
calculate_soft_scores <- function(spe, marker_stats, lineage_table,
                                  assay_name = "exprs",
                                  neg_strength = 0.8) {
  .assert_spe(spe)
  .assert_lineage_table(lineage_table)

  expr_mat <- SummarizedExperiment::assay(spe, assay_name)
  cell_types <- lineage_table$cell_type

  prob_mat <- matrix(0, nrow = ncol(spe), ncol = length(cell_types))
  colnames(prob_mat) <- cell_types
  rownames(prob_mat) <- colnames(spe)

  for (ct in cell_types) {
    pos <- lineage_table$pos_markers[lineage_table$cell_type == ct][[1]]
    neg <- lineage_table$neg_markers[lineage_table$cell_type == ct][[1]]

    pos <- intersect(pos, names(marker_stats))
    neg <- intersect(neg, names(marker_stats))
    if (length(pos) == 0) next

    score_sum <- numeric(ncol(spe))
    weight_sum <- 0

    for (m in pos) {
      st <- marker_stats[[m]]
      w <- if (is.finite(st$weight)) st$weight else 0.1
      x <- as.numeric(expr_mat[m, ])

      z <- (x - st$cutoff) / st$scale
      p <- stats::plogis(z)
      p[!is.finite(p)] <- 0

      score_sum  <- score_sum + p * w
      weight_sum <- weight_sum + w
    }

    base_score <- if (is.finite(weight_sum) && weight_sum > 0) score_sum / weight_sum else rep(0, ncol(spe))

    penalty <- rep(1, ncol(spe))
    if (length(neg) > 0) {
      for (m in neg) {
        st <- marker_stats[[m]]
        x <- as.numeric(expr_mat[m, ])
        z <- (x - st$cutoff) / st$scale
        p_bad <- stats::plogis(z)
        p_bad[!is.finite(p_bad)] <- 0
        penalty <- penalty * (1 - (p_bad * neg_strength))
      }
    }

    prob_mat[, ct] <- base_score * penalty
  }

  prob_mat
}
