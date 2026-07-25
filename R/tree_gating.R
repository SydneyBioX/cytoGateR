
#' Fit a 2-component Gaussian mixture model
#'
#' @param x numeric vector
#' @param cutoff_method Character. One of "mean" (midpoint of means) or
#'   "equal_posteriors" (threshold where component posteriors are equal).
#' @param gmm_model_names Optional character vector of model names to pass to
#'   [mclust::Mclust()] (e.g., "V" to forbid equal-variance in 1D).
#' @return list with GMM parameters or NULL
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, 0, 0.2), rnorm(50, 2, 0.2))
#' fit <- fit_gmm_2(x)
#' fit$cutoff
#' @export
fit_gmm_2 <- function(x,
                      cutoff_method = c("mean", "equal_posteriors"),
                      gmm_model_names = NULL) {
  cutoff_method <- match.arg(cutoff_method)
  x <- x[is.finite(x)]
  if (length(x) < 50L || length(unique(x)) < 3L) return(NULL)

  gmm <- tryCatch(
    mclust::Mclust(x, G = 2, modelNames = gmm_model_names, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(gmm) || length(gmm$parameters$mean) != 2) return(NULL)

  mu <- as.numeric(gmm$parameters$mean)
  # sig2 <- as.numeric(gmm$parameters$variance$sigmasq)

  # Check if sigmasq is a single value (Equal Variance model) and replicate it if so
  sig2 <- gmm$parameters$variance$sigmasq
  if (length(sig2) == 1) {
    sig2 <- rep(sig2, 2)
  }
  sig2 <- as.numeric(sig2)
  pi <- as.numeric(gmm$parameters$pro)

  ord <- order(mu)
  mu1 <- mu[ord[1]]; mu2 <- mu[ord[2]]
  s1 <- sqrt(sig2[ord[1]]); s2 <- sqrt(sig2[ord[2]])
  p1 <- pi[ord[1]]; p2 <- pi[ord[2]]

  sep <- abs(mu2 - mu1) / sqrt(s1^2 + s2^2)
  cutoff <- if (cutoff_method == "equal_posteriors") {
    gmm_equal_posterior_cutoff(mu1, mu2, s1, s2, p1, p2)
  } else {
    mean(c(mu1, mu2))
  }

  list(mu1 = mu1, mu2 = mu2, s1 = s1, s2 = s2, p1 = p1, p2 = p2,
       cutoff = cutoff, sep_score = sep)
}

#' Compute the cutoff where posteriors are equal
#'
#' @param mu1 Numeric mean of component 1.
#' @param mu2 Numeric mean of component 2.
#' @param s1 Numeric SD of component 1.
#' @param s2 Numeric SD of component 2.
#' @param p1 Numeric mixing weight of component 1.
#' @param p2 Numeric mixing weight of component 2.
#'
#' @return Numeric cutoff. Defaults to midpoint if root finding fails.
gmm_equal_posterior_cutoff <- function(mu1, mu2, s1, s2, p1, p2) {
  if (!all(is.finite(c(mu1, mu2, s1, s2, p1, p2))) || any(c(s1, s2) <= 0)) {
    return(mean(c(mu1, mu2)))
  }

  f <- function(x) {
    stats::dnorm(x, mean = mu1, sd = s1) * p1 -
      stats::dnorm(x, mean = mu2, sd = s2) * p2
  }

  lo <- min(mu1, mu2)
  hi <- max(mu1, mu2)
  root <- try(stats::uniroot(f, interval = c(lo, hi)), silent = TRUE)
  if (inherits(root, "try-error")) {
    return(mean(c(mu1, mu2)))
  }

  as.numeric(root$root)
}

#' Marker separability
#' @param x numeric
#' @return numeric
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, 0, 0.2), rnorm(50, 2, 0.2))
#' marker_separability(x)
#' @export
marker_separability <- function(x) {
  fit <- fit_gmm_2(x)
  if (is.null(fit)) return(NA_real_)
  fit$sep_score
}

#' Logistic marker score
#'
#' @param x numeric
#' @param cutoff numeric
#' @param scale numeric
#' @return numeric in [0,1]
#' @examples
#' score_marker_logistic(c(0, 1, 2), cutoff = 1, scale = 0.5)
#' @export
score_marker_logistic <- function(x, cutoff, scale) {
  stats::plogis((x - cutoff) / scale)
}

#' Rank based marker score
#'
#' @param x numeric
#' @return numeric in [0,1]
#' @examples
#' score_marker_rank(c(3, 1, 2, NA))
#' @export
score_marker_rank <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  (r - 1) / (sum(!is.na(r)) - 1 + 1e-9)
}

#' Build full-coverage gating tree
#'
#' @param expr_mat expression matrix
#' @param markers_pos positive markers
#' @param markers_neg negative markers
#' @param cell_idx indices
#' @param depth depth
#' @param max_depth max depth
#' @param min_cells minimum cells
#' @param min_score minimum separability
#' @param cutoff_method Cutoff selection method passed to [fit_gmm_2()].
#' @param gmm_model_names Optional character vector of model names to pass to
#'   [mclust::Mclust()] (e.g., "V" to forbid equal-variance in 1D).
#' @return tree object
#' @examples
#' set.seed(1)
#' expr_mat <- rbind(
#'   CD3 = c(rnorm(50, 0, 0.2), rnorm(50, 2, 0.2)),
#'   CD20 = c(rnorm(50, 2, 0.2), rnorm(50, 0, 0.2))
#' )
#' tree <- build_fullcoverage_tree(
#'   expr_mat,
#'   markers_pos = c("CD3", "CD20"),
#'   max_depth = 2,
#'   min_cells = 20,
#'   min_score = 0.1
#' )
#' tree$type
#' @export
build_fullcoverage_tree <- function(expr_mat,
                                    markers_pos,
                                    markers_neg = character(0),
                                    cell_idx = seq_len(ncol(expr_mat)),
                                    depth = 0,
                                    max_depth = 4,
                                    min_cells = 200,
                                    min_score = 0.5,
                                    cutoff_method = c("mean", "equal_posteriors"),
                                    gmm_model_names = NULL) {
  cutoff_method <- match.arg(cutoff_method)

  markers_pos <- intersect(markers_pos, rownames(expr_mat))
  markers_neg <- intersect(markers_neg, rownames(expr_mat))

  if (length(cell_idx) < min_cells || depth >= max_depth || length(markers_pos) == 0) {
    return(list(type = "leaf", depth = depth, cells = cell_idx))
  }

  Xpos <- as.matrix(expr_mat[markers_pos, cell_idx, drop = FALSE])
  pos_score <- colMeans(Xpos, na.rm = TRUE)

  neg_score <- 0
  if (length(markers_neg) > 0) {
    Xneg <- as.matrix(expr_mat[markers_neg, cell_idx, drop = FALSE])
    neg_score <- colMeans(Xneg, na.rm = TRUE)
  }

  target <- ifelse(pos_score - neg_score > median(pos_score - neg_score, na.rm = TRUE), 1, 0)

  marker_stats <- lapply(markers_pos, function(m) {
    x <- as.numeric(expr_mat[m, cell_idx])
    fit <- fit_gmm_2(
      x,
      cutoff_method = cutoff_method,
      gmm_model_names = gmm_model_names
    )
    if (is.null(fit)) return(c(marker = m, sep = NA, acc = NA, cutoff = NA, scale = NA))

    cutoff <- fit$cutoff
    scale <- sqrt(fit$s1^2 + fit$s2^2)
    if (!is.finite(scale) || scale == 0) scale <- .safe_mad(x)

    pred <- ifelse(x > cutoff, 1, 0)
    acc <- max(mean(pred == target, na.rm = TRUE),
               mean((1 - pred) == target, na.rm = TRUE))

    data.frame(marker = m, sep = fit$sep_score, acc = acc, cutoff = cutoff, scale = scale)
  })

  marker_stats <- do.call(rbind, marker_stats)
  marker_stats$sep <- as.numeric(marker_stats$sep)
  marker_stats$acc <- as.numeric(marker_stats$acc)
  marker_stats$combo <- marker_stats$sep + 2 * (marker_stats$acc - 0.5)

  marker_stats <- marker_stats[order(-marker_stats$combo), , drop = FALSE]
  best <- marker_stats[1, ]

  if (!is.finite(best$sep) || best$sep < min_score) {
    return(list(type = "leaf", depth = depth, cells = cell_idx))
  }

  best_marker <- best$marker
  cutoff <- as.numeric(best$cutoff)
  scale <- as.numeric(best$scale)

  x_best <- as.numeric(expr_mat[best_marker, cell_idx])
  left <- cell_idx[x_best <= cutoff]
  right <- cell_idx[x_best > cutoff]

  if (length(left) == 0 || length(right) == 0)
    return(list(type = "leaf", depth = depth, cells = cell_idx))

  remaining_pos <- setdiff(markers_pos, best_marker)

  list(
    type = "node",
    depth = depth,
    marker = best_marker,
    cutoff = cutoff,
    scale = scale,
    sep_score = as.numeric(best$sep),
    cells = cell_idx,
    left = build_fullcoverage_tree(expr_mat, remaining_pos, markers_neg, left,
                                   depth + 1, max_depth, min_cells, min_score,
                                   cutoff_method = cutoff_method,
                                   gmm_model_names = gmm_model_names),
    right = build_fullcoverage_tree(expr_mat, remaining_pos, markers_neg, right,
                                    depth + 1, max_depth, min_cells, min_score,
                                    cutoff_method = cutoff_method,
                                    gmm_model_names = gmm_model_names)
  )
}

#' Collect scores along a tree path
#'
#' Traverse a gating tree for a single cell and compute node-wise marker scores
#' along the deterministic path defined by comparing marker expression to each
#' node's cutoff.
#'
#' @param tree A gating tree as produced by [build_fullcoverage_tree()].
#' @param expr_mat Numeric matrix-like expression object with markers in rows and
#'   cells in columns (e.g. `assay(spe, "norm")`). Row names must include marker names.
#' @param cell_i Integer index of the cell/column to score.
#'
#' @return A numeric vector of per-node scores (may contain `NA_real_` if a marker value
#'   is not finite for that cell).
#' @examples
#' tree <- list(
#'   type = "node", marker = "CD3", cutoff = 1, scale = 0.5,
#'   left = list(type = "leaf"), right = list(type = "leaf")
#' )
#' expr_mat <- rbind(CD3 = c(0.5, 1.5))
#' collect_path_scores(tree, expr_mat, cell_i = 1)
#' @export
collect_path_scores <- function(tree, expr_mat, cell_i) {
  scores <- numeric(0)
  node <- tree

  while (!is.null(node) && node$type == "node") {
    m <- node$marker
    x <- as.numeric(expr_mat[m, cell_i])

    if (!is.finite(x)) {
      s <- NA_real_
    } else {
      s <- score_marker_logistic(x, node$cutoff, node$scale)
    }

    scores <- c(scores, s)
    if (is.finite(x) && x > node$cutoff) node <- node$right else node <- node$left
  }

  scores
}

#' Negative marker penalty
#'
#' Compute a multiplicative penalty factor in the range [0, 1] for a given cell
#' based on the expression of predefined negative markers.
#'
#' For each negative marker, a logistic scoring function is used to estimate
#' the probability (`p_bad`) that the marker is expressed above its expected
#' background level. The final penalty is computed as the product of
#' (1 - p_bad * neg_strength) across all negative markers.
#'
#' Higher negative-marker expression results in a stronger penalty (i.e.,
#' a smaller multiplicative factor).
#'
#' @param expr_mat A numeric matrix-like object containing expression values,
#'   with markers in rows and cells in columns.
#' @param neg_markers A character vector of marker names to penalize.
#' @param cell_i An integer index specifying the column (cell) to evaluate.
#' @param marker_stats A named list containing per-marker statistics. Each
#'   element must include at least:
#'   \describe{
#'     \item{cutoff}{Numeric value representing the logistic midpoint.}
#'     \item{scale}{Numeric value controlling the logistic slope.}
#'   }
#' @param neg_strength A numeric scalar in [0, 1] controlling the strength of
#'   the penalty. Larger values increase the impact of negative-marker expression.
#'
#' @return A numeric scalar in the range [0, 1] representing the multiplicative
#'   penalty factor for the specified cell.
#'
#' @details
#' If no valid negative markers are found in \code{marker_stats}, the function
#' returns 1 (no penalty).
#'
#' @examples
#' expr_mat <- rbind(CD3 = c(0.2, 2))
#' marker_stats <- list(CD3 = list(cutoff = 1, scale = 0.5))
#' neg_penalty(
#'   expr_mat,
#'   neg_markers = "CD3",
#'   cell_i = 1,
#'   marker_stats = marker_stats
#' )
#' @export
neg_penalty <- function(expr_mat, neg_markers, cell_i, marker_stats = NULL, neg_strength = 0.8) {
  neg_markers <- intersect(neg_markers, names(marker_stats))
  if (length(neg_markers) == 0) return(1)

  penalty <- 1
  for (m in neg_markers) {
    st <- marker_stats[[m]]
    x <- as.numeric(expr_mat[m, cell_i])

    if (is.finite(x)) {
      # Use existing logistic scoring function to get p_bad
      p_bad <- score_marker_logistic(x, st$cutoff, st$scale)
      penalty <- penalty * (1 - (p_bad * neg_strength))
    }
  }

  return(max(min(penalty, 1), 0))
}



#' Tree probability for one cell
#'
#' Compute the probability that a single cell belongs to a given cell type by
#' combining per-node scores along a tree path, applying a complexity-based
#' weight, and an optional negative marker penalty.
#'
#' @param tree A gating tree as produced by [build_fullcoverage_tree()].
#' @param expr_mat Numeric matrix-like expression object with markers in rows and
#'   cells in columns.
#' @param cell_i Integer index of the cell/column to score.
#' @param combine How to combine node-wise scores along the path. One of `"mean"`
#'   or `"product"`.
#' @param neg_markers Character vector of negative marker names used for penalty.
#' @param marker_stats A list of marker statistics (including cutoff and scale)
#'   produced by [fit_marker_stats()] used to calculate the penalty.
#' @param neg_strength Numeric scalar in [0,1] defining how strongly negative
#'   marker expression should penalize the final score. Default is 0.8.
#' @param lambda Numeric scalar representing the "Skepticism Factor." Higher
#'   values place a heavier penalty on cell types defined by fewer markers.
#'   Default is 0.2.
#'
#' @return Numeric scalar probability (typically in [0,1]) or `NA_real_` if no
#'   usable node scores are available for the cell.
#' @examples
#' tree <- list(
#'   type = "node", marker = "CD3", cutoff = 1, scale = 0.5,
#'   left = list(type = "leaf"), right = list(type = "leaf")
#' )
#' expr_mat <- rbind(CD3 = c(0.5, 1.5))
#' tree_prob(tree, expr_mat, cell_i = 2)
#' @export
tree_prob <- function(tree, expr_mat, cell_i,
                      combine = c("mean", "product"),
                      neg_markers = character(0),
                      marker_stats = NULL, # Add this
                      neg_strength = 0.8,
                      lambda = 2) { # Add this
  combine <- match.arg(combine)
  s <- collect_path_scores(tree, expr_mat, cell_i)
  s <- s[is.finite(s)]

  if (length(s) == 0) return(NA_real_)

  base <- if (combine == "mean") mean(s) else prod(s)

  n <- length(s)
  complexity_weight <- n / (n + lambda)
  # Pass marker_stats and strength to the new penalty function
  # base * neg_penalty(expr_mat, neg_markers, cell_i, marker_stats, neg_strength)
  # base * neg_penalty(expr_mat, neg_markers, cell_i)
  (base * complexity_weight) * neg_penalty(expr_mat, neg_markers, cell_i, marker_stats, neg_strength)
}
