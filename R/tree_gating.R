#' Fit a 2-component Gaussian mixture model
#'
#' @param x numeric vector
#' @return list with GMM parameters or NULL
#' @export
fit_gmm_2 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 50L || length(unique(x)) < 3L) return(NULL)

  gmm <- tryCatch(mclust::Mclust(x, G = 2, verbose = FALSE), error = function(e) NULL)
  if (is.null(gmm) || length(gmm$parameters$mean) != 2) return(NULL)

  mu <- as.numeric(gmm$parameters$mean)
  sig2 <- as.numeric(gmm$parameters$variance$sigmasq)
  pi <- as.numeric(gmm$parameters$pro)

  ord <- order(mu)
  mu1 <- mu[ord[1]]; mu2 <- mu[ord[2]]
  s1 <- sqrt(sig2[ord[1]]); s2 <- sqrt(sig2[ord[2]])
  p1 <- pi[ord[1]]; p2 <- pi[ord[2]]

  sep <- abs(mu2 - mu1) / sqrt(s1^2 + s2^2)
  cutoff <- mean(c(mu1, mu2))

  list(mu1 = mu1, mu2 = mu2, s1 = s1, s2 = s2, p1 = p1, p2 = p2,
       cutoff = cutoff, sep_score = sep)
}

#' Marker separability
#' @param x numeric
#' @return numeric
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
#' @export
score_marker_logistic <- function(x, cutoff, scale) {
  stats::plogis((x - cutoff) / scale)
}

#' Rank based marker score
#'
#' @param x numeric
#' @return numeric in [0,1]
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
#' @return tree object
#' @export
build_fullcoverage_tree <- function(expr_mat,
                                    markers_pos,
                                    markers_neg = character(0),
                                    cell_idx = seq_len(ncol(expr_mat)),
                                    depth = 0,
                                    max_depth = 4,
                                    min_cells = 200,
                                    min_score = 0.5) {

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
    fit <- fit_gmm_2(x)
    if (is.null(fit)) return(c(marker = m, sep = NA, acc = NA, cutoff = NA, scale = NA))

    cutoff <- fit$cutoff
    scale <- sqrt(fit$s1^2 + fit$s2^2)
    if (!is.finite(scale) || scale == 0) scale <- .safe_mad(x)

    pred <- ifelse(x > cutoff, 1, 0)
    acc <- max(mean(pred == target, na.rm = TRUE),
               mean((1 - pred) == target, na.rm = TRUE))

    c(marker = m, sep = fit$sep_score, acc = acc, cutoff = cutoff, scale = scale)
  })

  marker_stats <- as.data.frame(do.call(rbind, marker_stats))
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
                                   depth + 1, max_depth, min_cells, min_score),
    right = build_fullcoverage_tree(expr_mat, remaining_pos, markers_neg, right,
                                    depth + 1, max_depth, min_cells, min_score)
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
#' Compute a multiplicative penalty in [0,1] based on the mean expression of
#' negative markers for a given cell. Higher negative-marker expression leads
#' to stronger penalty (smaller value).
#'
#' @param expr_mat Numeric matrix-like expression object with markers in rows and
#'   cells in columns.
#' @param neg_markers Character vector of marker names to penalize.
#' @param cell_i Integer index of the cell/column to evaluate.
#'
#' @return A numeric scalar in [0,1] representing the penalty factor.
#' @export
neg_penalty <- function(expr_mat, neg_markers, cell_i) {
  neg_markers <- intersect(neg_markers, rownames(expr_mat))
  if (length(neg_markers) == 0) return(1)

  x <- as.numeric(expr_mat[neg_markers, cell_i])
  x <- x[is.finite(x)]
  if (length(x) == 0) return(1)

  p <- 1 - mean(x, na.rm = TRUE)
  max(min(p, 1), 0)
}

#' Tree probability for one cell
#'
#' Compute the probability that a single cell belongs to a given cell type by
#' combining per-node scores along a tree path and applying an optional negative
#' marker penalty.
#'
#' @param tree A gating tree as produced by [build_fullcoverage_tree()].
#' @param expr_mat Numeric matrix-like expression object with markers in rows and
#'   cells in columns.
#' @param cell_i Integer index of the cell/column to score.
#' @param combine How to combine node-wise scores along the path. One of `"mean"`
#'   or `"product"`.
#' @param neg_markers Character vector of negative marker names used for penalty.
#'
#' @return Numeric scalar probability (typically in [0,1]) or `NA_real_` if no
#'   usable node scores are available for the cell.
#' @export
tree_prob <- function(tree, expr_mat, cell_i,
                      combine = c("mean", "product"),
                      neg_markers = character(0)) {
  combine <- match.arg(combine)
  s <- collect_path_scores(tree, expr_mat, cell_i)
  s <- s[is.finite(s)]

  if (length(s) == 0) return(NA_real_)

  base <- if (combine == "mean") mean(s) else prod(s)
  base * neg_penalty(expr_mat, neg_markers, cell_i)
}

