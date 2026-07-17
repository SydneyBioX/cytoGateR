#' Auto-detect distribution shape and select best fitting method
#'
#' @param x Numeric vector of expression values.
#' @param zero_thresh Value below which a cell is considered "near zero".
#' @param zero_inflate_frac If more than this fraction are near zero, flag as zero-inflated.
#' @param skew_thresh If skewness exceeds this, flag as heavily skewed.
#' @return Character string: one of "gmm_log1p", "gmm", "gamma", "norm"
detect_fit_method <- function(x,
                              zero_thresh       = 0.1,
                              zero_inflate_frac = 0.4,
                              skew_thresh       = 2.0) {
  x <- x[is.finite(x) & x >= 0]

  frac_zero <- mean(x < zero_thresh)
  skewness  <- mean((x - mean(x))^3) / (sd(x)^3 + 1e-10)
  all_pos   <- all(x > 0)

  if (frac_zero > zero_inflate_frac || abs(skewness) > skew_thresh) {
    # Heavy zero-inflation or strong skew → log1p then Gaussian (mclust)
    return("gmm_log1p")
  } else if (all_pos && abs(skewness) > 1.0) {
    # Moderately skewed, strictly positive → gamma worth trying
    return("gamma")
  } else if (abs(skewness) <= 1.0) {
    # Roughly symmetric → standard Gaussian mixture (mclust)
    return("gmm")
  } else {
    # Fallback: symmetric but not strictly positive → mixtools normal
    return("norm")
  }
}


#' Fit best mixture model with automatic method detection
#'
#' @param x Numeric vector.
#' @param cutoff_method One of "mean" or "equal_posteriors".
#' @param gmm_model_names Optional mclust model names.
#' @param force_method If not NULL, skip detection and use this method directly.
#'   One of "gmm", "gmm_log1p", "gamma", "norm".
#' @param verbose If TRUE, message which method was selected.
#' @return list with mixture parameters (same fields as fit_gmm_2) or NULL
fit_auto_2 <- function(x,
                       cutoff_method   = c("mean", "equal_posteriors"),
                       gmm_model_names = NULL,
                       force_method    = NULL,
                       verbose         = TRUE) {
  cutoff_method <- match.arg(cutoff_method)
  x_clean <- x[is.finite(x)]

  # --- Method Selection ---
  method <- if (!is.null(force_method)) {
    force_method
  } else {
    detect_fit_method(x_clean)
  }

  if (verbose) message(sprintf("  [fit_auto_2] selected method: %s", method))

  # --- Helper: back-transform a log1p fit to original scale ---
  .backtransform_log1p_fit <- function(fit_log) {
    if (is.null(fit_log)) return(NULL)
    mu1_log <- fit_log$mu1
    mu2_log <- fit_log$mu2
    s1_log  <- fit_log$s1
    s2_log  <- fit_log$s2
    fit_log$cutoff <- expm1(fit_log$cutoff)
    fit_log$mu1    <- expm1(mu1_log)
    fit_log$mu2    <- expm1(mu2_log)
    fit_log$s1     <- exp(mu1_log) * s1_log  # delta method
    fit_log$s2     <- exp(mu2_log) * s2_log  # delta method
    fit_log
  }

  # --- Helper: fit gmm on log1p scale and back-transform ---
  .fit_gmm_log1p <- function() {
    fit_log <- fit_gmm_2(log1p(x_clean),
                         cutoff_method   = cutoff_method,
                         gmm_model_names = gmm_model_names)
    .backtransform_log1p_fit(fit_log)
  }

  # --- Helper: fit norm on log1p scale and back-transform ---
  .fit_norm_log1p <- function() {
    fit_log <- fit_norm_mix_2(log1p(x_clean),
                              cutoff_method = cutoff_method)
    .backtransform_log1p_fit(fit_log)
  }

  # --- Fit ---
  fit <- switch(method,

                "gmm" = fit_gmm_2(x_clean,
                                  cutoff_method   = cutoff_method,
                                  gmm_model_names = gmm_model_names),

                "gmm_log1p" = .fit_gmm_log1p(),

                "gamma" = {
                  fit_g <- fit_gamma_mix_2(x_clean, cutoff_method = cutoff_method)
                  if (is.null(fit_g)) {
                    if (verbose) message("  [fit_auto_2] gamma diverged, falling back to gmm_log1p")
                    .fit_gmm_log1p()
                  } else {
                    fit_g
                  }
                },

                "norm" = fit_norm_mix_2(x_clean, cutoff_method = cutoff_method),

                "norm_log1p" = .fit_norm_log1p(),

                stop(sprintf("Unknown method '%s'. Choose one of: 'gmm', 'gmm_log1p', 'gamma', 'norm', 'norm_log1p'.", method))
  )

  return(fit)
}

#' Fit a 2-component Gamma mixture model
#'
#' @param x numeric vector
#' @param cutoff_method Character. One of "mean" (midpoint of means) or
#'   "equal_posteriors" (threshold where component posteriors are equal).
#' @param maxit Integer. Maximum iterations for EM algorithm. Default 1000.
#' @param epsilon Numeric. Convergence tolerance. Default 1e-4.
#' @return list with GMM parameters or NULL
fit_gamma_mix_2 <- function(x,
                            cutoff_method = c("mean", "equal_posteriors"),
                            maxit = 1000,
                            epsilon = 1e-4) {
  cutoff_method <- match.arg(cutoff_method)
  x <- x[is.finite(x) & x > 0]  # gamma requires strictly positive values
  if (length(x) < 50L || length(unique(x)) < 3L) return(NULL)

  fit <- tryCatch(
    mixtools::gammamixEM(x, k = 2, maxit = maxit, epsilon = epsilon),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  # gamma.pars: matrix with rows = shape (alpha) and rate (beta), cols = components
  shapes <- fit$gamma.pars[1, ]  # alpha
  rates  <- fit$gamma.pars[2, ]  # beta
  pi     <- fit$lambda            # mixing proportions

  # Component means: mean of gamma = shape / rate
  means <- shapes / rates

  ord <- order(means)
  shape1 <- shapes[ord[1]]; shape2 <- shapes[ord[2]]
  rate1  <- rates[ord[1]];  rate2  <- rates[ord[2]]
  p1     <- pi[ord[1]];     p2     <- pi[ord[2]]

  # Mean and SD of each gamma component
  mu1 <- shape1 / rate1;  mu2 <- shape2 / rate2
  s1  <- sqrt(shape1) / rate1;  s2 <- sqrt(shape2) / rate2

  sep <- abs(mu2 - mu1) / sqrt(s1^2 + s2^2)

  cutoff <- if (cutoff_method == "equal_posteriors") {
    gamma_mix_equal_posterior_cutoff(shape1, rate1, shape2, rate2, p1, p2, mu1, mu2)
  } else {
    mean(c(mu1, mu2))
  }

  list(mu1 = mu1, mu2 = mu2, s1 = s1, s2 = s2,
       shape1 = shape1, rate1 = rate1,
       shape2 = shape2, rate2 = rate2,
       p1 = p1, p2 = p2,
       cutoff = cutoff, sep_score = sep)
}

#' Compute the cutoff where gamma mixture posteriors are equal
#'
#' @param shape1 Numeric shape of component 1.
#' @param rate1  Numeric rate of component 1.
#' @param shape2 Numeric shape of component 2.
#' @param rate2  Numeric rate of component 2.
#' @param p1 Numeric mixing weight of component 1.
#' @param p2 Numeric mixing weight of component 2.
#' @param mu1 Numeric mean of component 1 (used as fallback midpoint).
#' @param mu2 Numeric mean of component 2 (used as fallback midpoint).
#'
#' @return Numeric cutoff. Defaults to midpoint of means if root finding fails.
gamma_mix_equal_posterior_cutoff <- function(shape1, rate1, shape2, rate2,
                                             p1, p2, mu1, mu2) {
  if (!all(is.finite(c(shape1, rate1, shape2, rate2, p1, p2))) ||
      any(c(shape1, rate1, shape2, rate2) <= 0)) {
    return(mean(c(mu1, mu2)))
  }
  f <- function(x) {
    stats::dgamma(x, shape = shape1, rate = rate1) * p1 -
      stats::dgamma(x, shape = shape2, rate = rate2) * p2
  }
  lo <- min(mu1, mu2)
  hi <- max(mu1, mu2)
  root <- try(stats::uniroot(f, interval = c(lo, hi)), silent = TRUE)
  if (inherits(root, "try-error")) return(mean(c(mu1, mu2)))
  as.numeric(root$root)
}


#' Fit a 2-component Normal mixture model via mixtools
#'
#' @param x numeric vector
#' @param cutoff_method Character. One of "mean" (midpoint of means) or
#'   "equal_posteriors" (threshold where component posteriors are equal).
#' @param maxit Integer. Maximum iterations for EM algorithm. Default 1000.
#' @param epsilon Numeric. Convergence tolerance. Default 1e-4.
#' @return list with mixture parameters or NULL
fit_norm_mix_2 <- function(x,
                           cutoff_method = c("mean", "equal_posteriors"),
                           maxit = 1000,
                           epsilon = 1e-4) {
  cutoff_method <- match.arg(cutoff_method)
  x <- x[is.finite(x)]
  if (length(x) < 50L || length(unique(x)) < 3L) return(NULL)

  fit <- tryCatch(
    mixtools::normalmixEM(x, k = 2, maxit = maxit, epsilon = epsilon),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  mu  <- fit$mu     # means
  sig <- fit$sigma  # SDs
  pi  <- fit$lambda # mixing proportions

  ord <- order(mu)
  mu1 <- mu[ord[1]]; mu2 <- mu[ord[2]]
  s1  <- sig[ord[1]]; s2 <- sig[ord[2]]
  p1  <- pi[ord[1]];  p2 <- pi[ord[2]]

  sep <- abs(mu2 - mu1) / sqrt(s1^2 + s2^2)

  cutoff <- if (cutoff_method == "equal_posteriors") {
    gmm_equal_posterior_cutoff(mu1, mu2, s1, s2, p1, p2)  # reuse your existing function
  } else {
    mean(c(mu1, mu2))
  }

  list(mu1 = mu1, mu2 = mu2, s1 = s1, s2 = s2,
       p1 = p1, p2 = p2,
       cutoff = cutoff, sep_score = sep)
}


#' Marker separability
#' @param x numeric
#' @return numeric
#' @export

marker_separability <- function(x,
                                 fit_method = "em_norm",
                                 cutoff_method = c("mean", "equal_posteriors"),
                                 gmm_model_names = NULL) {
  fit_method <- match.arg(fit_method)

  if(fit_method=="auto"){
    fit <- fit_auto_2(
      x,
      cutoff_method   = cutoff_method,
      gmm_model_names = gmm_model_names,
      verbose         = TRUE           # set FALSE to silence per-marker messages
    )

  }else{
    fit <- switch(fit_method,
                  "em_mclust" = fit_gmm_2(x,
                                          cutoff_method = cutoff_method,
                                          gmm_model_names = gmm_model_names),
                  "em_norm"   = fit_norm_mix_2(x,
                                               cutoff_method = cutoff_method),
                  # "em_gamma"  = fit_gamma_mix_2(x,
                  #                            cutoff_method = cutoff_method)
    )
  }

  if (is.null(fit)) return(NA_real_)
  fit$sep_score
}




#' Multi-Metric Weighted kNN
#'
#' @param train_data Data frame or matrix of training cells (rows = markers, cols = cells).
#' @param test_data Data frame or matrix of test cells (rows = markers, cols = cells).
#' @param train_labels Factor of labels for the training data.
#' @param k Number of neighbors.
#' @param method One of "pearson", "spearman", "cosine", or "euclidean".
#'
#' @return A list containing predicted 'labels' and 'probs'.
predict_wknn_multi <- function(train_data,
                               test_data,
                               train_labels,
                               k = 5,
                               method = "pearson",
                               return_matrix = FALSE) {

  train_mat <- t(as.matrix(train_data))
  test_mat  <- t(as.matrix(test_data))

  if (nrow(train_mat) != nrow(test_mat)) {
    stop(sprintf(
      "Dimension mismatch! Train has %d markers, Test has %d markers.",
      nrow(train_mat), nrow(test_mat)
    ))
  }

  all_classes <- levels(train_labels)

  # 1. Similarity/Distance Matrix
  if (method %in% c("pearson", "spearman")) {
    score_mat  <- cor(test_mat, train_mat, method = method)
    # cor() can return NA if a cell has zero variance — replace with 0
    score_mat[is.na(score_mat)] <- 0
    is_distance <- FALSE
  } else if (method == "cosine") {
    cp         <- crossprod(test_mat, train_mat)
    rn         <- sqrt(colSums(test_mat^2))
    cn         <- sqrt(colSums(train_mat^2))
    score_mat  <- cp / outer(rn, cn)
    score_mat[is.na(score_mat)] <- 0
    is_distance <- FALSE
  } else if (method == "euclidean") {
    score_mat  <- as.matrix(proxy::dist(t(test_mat), t(train_mat), method = "Euclidean"))
    is_distance <- TRUE
  } else {
    stop(sprintf("Unknown method: '%s'. Use pearson, spearman, cosine, or euclidean.", method))
  }

  k <- min(k, ncol(train_mat))  # guard against k > n_train_cells

  # 2. Per-cell neighbour voting
  results <- apply(score_mat, 1, function(scores) {

    if (is_distance) {
      top_k_idx <- order(scores, decreasing = FALSE)[seq_len(k)]
      weights   <- 1 / (scores[top_k_idx] + 1e-6)
    } else {
      top_k_idx <- order(scores, decreasing = TRUE)[seq_len(k)]
      weights   <- pmax(scores[top_k_idx], 0)
    }

    top_labels  <- train_labels[top_k_idx]
    total_weight <- sum(weights)

    # Initialise probability vector over ALL known classes
    prob_dist <- setNames(numeric(length(all_classes)), all_classes)

    if (total_weight > 0) {
      label_sums <- tapply(weights, top_labels, sum)
      prob_dist[names(label_sums)] <- label_sums / total_weight
      best_label <- names(which.max(label_sums))   # which.max is safer than sort()[1]
      best_prob  <- max(label_sums) / total_weight
    } else {
      # All weights zero: fall back to plurality of hard votes
      vote_table <- table(top_labels)
      best_label <- names(which.max(vote_table))
      best_prob  <- 0
    }

    list(label = best_label, prob = best_prob, dist = prob_dist)
  })

  # 3. Collate output
  out <- list(
    labels = sapply(results, `[[`, "label"),
    probs  = sapply(results, `[[`, "prob")
  )

  if (return_matrix) {
    out$prob_matrix <- do.call(rbind, lapply(results, `[[`, "dist"))
    rownames(out$prob_matrix) <- rownames(score_mat)
  }

  return(out)
}


#' Recursive Hierarchical kNN with BiocParallel Ensemble
#'
#' @param spe SpatialExperiment object.
#' @param hier_ref List of node-specific references from build_hierarchical_reference.
#' @param hc_tree hclust object defining the hierarchy.
#' @param assay_name Assay to use (default "exprs").
#' @param threshold Confidence threshold (average kNN probability).
#' @param agreement_threshold Minimum \% of ensemble members that must agree (0-1).
#' @param k Number of neighbors for kNN.
#' @param repeats Number of bootstrap repeats per distance method.
#' @param dist_methods Vector of methods, e.g., c("pearson", "cosine").
#' @param BPPARAM BiocParallel parameter (default: SerialParam()).
#' @param out_col Column name for final labels.
#' @export
predict_hierarchical_knn_recursive <- function(spe,
                                               hier_ref,
                                               hc_tree,
                                               assay_name = "exprs",
                                               threshold = 0.7,
                                               agreement_threshold = 0.8,
                                               k = 5,
                                               repeats = 5,
                                               dist_methods = c("pearson", "cosine"),
                                               BPPARAM = BiocParallel::SerialParam(),
                                               out_col = "hier_label") {

  .assert_spe(spe)
  n_cells <- ncol(spe)
  feat_mat <- SummarizedExperiment::assay(spe, assay_name)

  # Initialize results in parent scope
  final_labels <- rep(NA_character_, n_cells)
  root_node_idx <- nrow(hc_tree$merge)

  # --- Internal Recursive Processor ---
  process_node <- function(node_idx, active_indices) {
    if (length(active_indices) == 0) return()

    node_id <- paste0("Node_", node_idx)
    ref <- hier_ref[[node_id]]

    # Subset features for this node (Cells x Markers)
    test_data <- t(as.matrix(feat_mat[ref$markers, active_indices, drop = FALSE]))
    train_data <- as.matrix(ref$train_data) # Already Cells x Markers

    # --- Ensemble Step ---
    # Run multiple methods and bootstraps in parallel
    task_grid <- expand.grid(method = dist_methods, r = seq_len(repeats), stringsAsFactors = FALSE)

    ensemble_results <- BiocParallel::bplapply(seq_len(nrow(task_grid)), function(i) {
      m <- task_grid$method[i]
      # Bootstrap 80% of training data
      boot_idx <- sample(seq_len(nrow(train_data)), size = floor(0.8 * nrow(train_data)))

      # Use our multi-metric engine
      predict_wknn_multi(
        train_data = train_data[boot_idx, , drop = FALSE],
        test_data = test_data,
        train_labels = ref$train_labels[boot_idx],
        k = k,
        method = m
      )
    }, BPPARAM = BPPARAM)

    # --- Consensus Gathering ---
    all_votes <- do.call(cbind, lapply(ensemble_results, `[[`, "labels"))
    all_probs <- do.call(cbind, lapply(ensemble_results, `[[`, "probs"))

    # 1. Majority Vote
    node_preds <- apply(all_votes, 1, function(x) {
      tbl <- table(x)
      names(sort(tbl, decreasing = TRUE))[1]
    })

    # 2. Agreement Score & Average Confidence
    node_agreement <- rowSums(all_votes == node_preds) / ncol(all_votes)
    node_avg_probs <- rowMeans(all_probs)

    # --- Gatekeeping ---
    # Cell must pass BOTH the probability threshold and the ensemble agreement
    uncertain_mask <- (node_avg_probs < threshold) | (node_agreement < agreement_threshold)

    if (any(uncertain_mask)) {
      final_labels[active_indices[uncertain_mask]] <<- paste0(node_id, "_unassigned")
    }

    # --- Routing ---
    certain_idx <- which(!uncertain_mask)
    if (length(certain_idx) > 0) {
      preds_certain <- node_preds[certain_idx]
      indices_certain <- active_indices[certain_idx]

      for (choice in c("Left", "Right")) {
        group_idx <- indices_certain[preds_certain == choice]
        if (length(group_idx) == 0) next

        side_idx <- if(choice == "Left") 1 else 2
        child_val <- hc_tree$merge[node_idx, side_idx]

        if (child_val < 0) {
          final_labels[group_idx] <<- hc_tree$labels[-child_val]
        } else {
          process_node(child_val, group_idx)
        }
      }
    }
  }

  message(sprintf("Starting recursive ensemble classification for %d cells...", n_cells))
  process_node(root_node_idx, seq_len(n_cells))

  SummarizedExperiment::colData(spe)[[out_col]] <- final_labels
  return(spe)
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

    c(marker = m, sep = fit$sep_score, acc = acc, cutoff = cutoff, scale = scale)
  })


  # marker_stats <- lapply(markers_pos, function(m) {
  #   x <- as.numeric(expr_mat[m, cell_idx])
  #
  #   fit <- switch(fit_method,
  #                 "em_mclust"   = fit_gmm_2(x,
  #                                     cutoff_method = cutoff_method,
  #                                     gmm_model_names = gmm_model_names),
  #                 # "em_gamma" = fit_gamma_mix_2(x,
  #                 #                           cutoff_method = cutoff_method),
  #                 "em_norm"  = fit_norm_mix_2(x,
  #                                          cutoff_method = cutoff_method),
  #                 stop("fit_method must be one of 'em_mclust', or 'em_norm'")
  #   )
  #
  #   if (is.null(fit)) return(c(marker = m, sep = NA, acc = NA, cutoff = NA, scale = NA))
  #
  #   cutoff <- fit$cutoff
  #   scale  <- sqrt(fit$s1^2 + fit$s2^2)  # works unchanged — s1/s2 are SDs in all three
  #   if (!is.finite(scale) || scale == 0) scale <- .safe_mad(x)
  #
  #   pred <- ifelse(x > cutoff, 1, 0)
  #   acc  <- max(mean(pred == target, na.rm = TRUE),
  #               mean((1 - pred) == target, na.rm = TRUE))
  #
  #   c(marker = m, sep = fit$sep_score, acc = acc, cutoff = cutoff, scale = scale)
  # })

  #
  # marker_stats <- lapply(markers_pos, function(m) {
  #   x <- as.numeric(expr_mat[m, cell_idx])

  # fit <- if (fit_method == "auto") {
  #   fit_auto_2(x,
  #              cutoff_method   = cutoff_method,
  #              gmm_model_names = gmm_model_names,
  #              verbose         = FALSE)
  # } else {
  # fit <- switch(fit_method,
  #        "em_mclust" = fit_gmm_2(x,
  #                                cutoff_method   = cutoff_method,
  #                                gmm_model_names = gmm_model_names),
  #        "em_norm"   = fit_norm_mix_2(x,
  #                                     cutoff_method = cutoff_method),
  #        "em_gamma"  = fit_gamma_mix_2(x,
  #                                      cutoff_method = cutoff_method),
  #        stop("fit_method must be one of 'auto', 'em_mclust', 'em_norm', 'em_gamma'")
  # )
  # }
  #
  #   if (is.null(fit)) return(c(marker = m, sep = NA, acc = NA, cutoff = NA, scale = NA))
  #
  #   cutoff <- fit$cutoff
  #   scale  <- sqrt(fit$s1^2 + fit$s2^2)
  #   if (!is.finite(scale) || scale == 0) scale <- .safe_mad(x)
  #
  #   pred <- ifelse(x > cutoff, 1, 0)
  #   acc  <- max(mean(pred == target, na.rm = TRUE),
  #               mean((1 - pred) == target, na.rm = TRUE))
  #
  #   c(marker = m, sep = fit$sep_score, acc = acc, cutoff = cutoff, scale = scale)
  # })

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
                                   depth + 1, max_depth, min_cells, min_score,
                                   cutoff_method = cutoff_method,
                                   gmm_model_names = gmm_model_names),
    right = build_fullcoverage_tree(expr_mat, remaining_pos, markers_neg, right,
                                    depth + 1, max_depth, min_cells, min_score,
                                    cutoff_method = cutoff_method,
                                    gmm_model_names = gmm_model_names)
  )
}


#' Train kNN Reference for Cell Type Prediction
#'
#' Trains a kNN-based cell type classifier using a consensus cleaning approach.
#' Core cells are iteratively validated across repeated cross-validation folds,
#' and cells with low label agreement are reassigned as unknown before a clean
#' reference set is returned for downstream prediction.
#'
#' @param spe A \code{SpatialExperiment} or \code{SummarizedExperiment} object
#'   containing cell data.
#' @param label_col Character string specifying the column in \code{colData(spe)}
#'   containing the initial cell type labels. Default is \code{"cutoff_label"}.
#' @param assay_name Character string specifying the assay to use for marker
#'   expression. Default is \code{"norm"}.
#' @param unknown_label Character string specifying the label used to identify
#'   unknown or unclassified cells, which are excluded from training. Default is
#'   \code{"Unknown"}.
#' @param features Either \code{"all"} to use all features in the assay, or a
#'   character vector of feature names to subset. Default is \code{"all"}.
#' @param cv_folds Integer specifying the number of cross-validation folds used
#'   during consensus cleaning. Default is \code{5}.
#' @param repeats Integer specifying the number of cleaning repetitions. Higher
#'   values produce more stable agreement rates. Default is \code{10}.
#' @param agreement_thresh Numeric value between 0 and 1 specifying the minimum
#'   proportion of repeats a cell must be correctly predicted to be retained in
#'   the reference set. Default is \code{0.8}.
#' @param seed Optional integer for random seed to ensure reproducibility.
#'   Default is \code{NULL}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{spe}}{The input \code{SpatialExperiment} object with a new
#'     \code{cleaned_core_label} column added to \code{colData}.}
#'   \item{\code{reference_data}}{A data frame of marker expression for the
#'     high-confidence reference cells.}
#'   \item{\code{reference_labels}}{A factor of cell type labels corresponding
#'     to \code{reference_data}.}
#'   \item{\code{features_used}}{Character vector of feature names used in
#'     training.}
#'   \item{\code{agreement_rates}}{Named numeric vector of per-cell agreement
#'     rates from consensus cleaning.}
#'   \item{\code{cleaned_core_names}}{Character vector of cell names retained
#'     in the reference set after cleaning.}
#' }
#'
#' @export
train_custom_knn <- function(spe,
                             label_col = "cutoff_label",
                             assay_name = "norm",
                             unknown_label = "Unknown",
                             features = "all",
                             cv_folds = 5,
                             repeats = 10,
                             agreement_thresh = 0.8,
                             seed = NULL) {

  if (!requireNamespace("class", quietly = TRUE)) stop("Please install 'class' package.")
  if (!is.null(seed)) set.seed(seed)
  .assert_spe(spe)

  # 1. Feature Prep
  feat_mat <- SummarizedExperiment::assay(spe, assay_name)
  features_use <- if (length(features) == 1 && features == "all") rownames(feat_mat) else intersect(features, rownames(feat_mat))

  lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
  feature_df <- as.data.frame(t(feat_mat[features_use, , drop = FALSE]))

  core_idx <- which(lab_vec != unknown_label & !is.na(lab_vec))
  if (length(core_idx) < 10) stop("Not enough core cells for CV cleaning.")

  core_df <- feature_df[core_idx, ]
  core_labels <- factor(lab_vec[core_idx])

  # --- STAGE 1: CONSENSUS CLEANING ---
  message(sprintf("Starting kNN label cleaning: %d repeats...", repeats))
  match_counts <- setNames(numeric(nrow(core_df)), rownames(core_df))

  for (r in seq_len(repeats)) {
    # Stratified fold assignment
    fold_assign <- integer(nrow(core_df))
    for (cls in levels(core_labels)) {
      cls_idx <- which(core_labels == cls)
      fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
    }

    for (k in seq_len(cv_folds)) {
      train_idx <- which(fold_assign != k); test_idx <- which(fold_assign == k)

      preds <- class::knn(train = core_df[train_idx, ],
                          test = core_df[test_idx, ],
                          cl = core_labels[train_idx], k = 5)

      match_counts[test_idx] <- match_counts[test_idx] + as.numeric(preds == core_labels[test_idx])
    }
  }

  agreement_rate <- match_counts / repeats

  # 2. Update labels in the SPE object (Matches your RF Logic)
  inconsistent_names <- names(agreement_rate)[agreement_rate < agreement_thresh]
  valid_core_names <- names(agreement_rate)[agreement_rate >= agreement_thresh]

  cleaned_labels <- lab_vec
  # Match by cell names to ensure accuracy
  cleaned_labels[colnames(spe) %in% inconsistent_names] <- unknown_label
  SummarizedExperiment::colData(spe)$cleaned_core_label <- cleaned_labels

  message(sprintf("Cleaning complete. Removed %d inconsistent core cells.", length(inconsistent_names)))

  # 3. Return List (Structure matches train_custom_randomForest)
  return(list(
    spe = spe,                                      # Now includes the new column
    reference_data = core_df[valid_core_names, ],    # The "Golden Set" for prediction
    reference_labels = core_labels[rownames(core_df) %in% valid_core_names],
    features_used = features_use,
    agreement_rates = agreement_rate,
    cleaned_core_names = valid_core_names
  ))
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




#' Train a custom random forest with consensus label cleaning
#'
#' Train a \code{ranger} random forest classifier on labelled cells from a
#' \code{SpatialExperiment} or \code{SingleCellExperiment}. The workflow has two
#' stages:
#' \enumerate{
#'   \item \strong{Consensus cleaning:} repeated stratified cross-validation is
#'   used to estimate how consistently each labelled cell is predicted as its
#'   original label. Cells with agreement below \code{agreement_thresh} are
#'   flagged as unreliable and relabeled as \code{unknown_label} in
#'   \code{colData(spe)$cleaned_core_label}.
#'   \item \strong{Model evaluation and final training:} cross-validation is run
#'   on the cleaned labelled set to report accuracy and macro-F1, and a final
#'   probability random forest model is trained on all cleaned labelled cells.
#' }
#'
#' @param spe A \code{SpatialExperiment} or \code{SingleCellExperiment}.
#' @param label_col Character scalar naming the label column in
#'   \code{colData(spe)} (default \code{"cutoff_label"}).
#' @param assay_name Character scalar naming the assay used as features
#'   (default \code{"norm"}).
#' @param unknown_label Character label treated as unlabeled and excluded from
#'   training/cleaning (default \code{"Unknown"}).
#' @param num.trees Integer number of trees for the final \code{ranger} model
#'   (default \code{200}).
#' @param mtry Optional integer \code{mtry} for \code{ranger}. If \code{NULL},
#'   uses \code{floor(sqrt(p))} where \code{p} is the number of features used.
#' @param seed Optional integer seed for reproducibility.
#' @param features Character vector of feature (marker) names to use, or
#'   \code{"all"} to use all assay rows (default \code{"all"}).
#' @param cv_folds Integer number of folds used for stratified cross-validation
#'   (default \code{5}).
#' @param repeats Integer number of repeated CV rounds used during consensus
#'   cleaning (default \code{10}).
#' @param agreement_thresh Numeric in [0,1] specifying the minimum agreement rate
#'   required to keep an originally labelled cell during cleaning (default
#'   \code{0.8}).
#' @param num_threads Integer requested number of threads for model fitting.
#'
#' @return A named list with components:
#' \describe{
#'   \item{spe}{The input \code{spe} with an added \code{colData} column
#'     \code{cleaned_core_label}, where low-agreement labelled cells are set to
#'     \code{unknown_label}.}
#'   \item{model}{A fitted \code{ranger} model trained on the cleaned labelled
#'     cells with \code{probability = TRUE}.}
#'   \item{agreement_rates}{Named numeric vector of per-cell agreement rates from
#'     the consensus cleaning stage (indexed by training-row names).}
#'   \item{features_used}{Character vector of feature names used for training.}
#'   \item{metrics}{List of evaluation outputs on the cleaned labelled set:
#'     \describe{
#'       \item{confusion_matrix}{Confusion matrix for pooled CV predictions.}
#'       \item{accuracy}{Overall pooled CV accuracy.}
#'       \item{f1_macro}{Overall pooled macro-F1.}
#'       \item{cv_overall}{Per-fold data.frame with accuracy and macro-F1.}
#'       \item{cv_class_metrics}{Per-fold, per-class precision/recall/F1 table.}
#'     }}
#'   \item{test_pred}{Factor of pooled cross-validated predictions on cleaned labelled cells.}
#'   \item{test_truth}{Factor of pooled cross-validated true labels on cleaned labelled cells.}
#' }
#'
#' @details
#' During consensus cleaning, models are trained without probability estimation.
#' The final model is trained with \code{probability = TRUE} to support downstream
#' thresholding workflows.
#'
#' @export
# train_custom_rf <- function(spe,
#                             label_col = "custom_label",
#                             assay_name = "norm",
#                             unknown_label = "Unknown",
#                             train_frac = 0.8,
#                             num.trees = 200,
#                             mtry = NULL,
#                             seed = NULL,
#                             features="lineage",
#                             cv_folds = 5) {
#   .assert_spe(spe)
#   if (!label_col %in% colnames(SummarizedExperiment::colData(spe))) {
#     stop("label_col not found in colData(spe).")
#   }
#   if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
#     stop(paste0("Assay '", assay_name, "' not found in spe."))
#   }
#
#
#
#   if (is.character(features)) {
#
#     if (length(features) == 1 && features %in% c("all")) {
#       # keyword mode
#       features_use <- rownames(feat_mat)
#     } else {
#       # explicit marker vector
#       features_use <- features
#     }
#
#   } else {
#     stop("features must be 'lineage', 'all', or a character vector of markers")
#   }
#
#
#
#
#   lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
#   feat_mat <- SummarizedExperiment::assay(spe, assay_name)
#
#
#
#
#   features_use <- unique(features_use)
#   present <- intersect(features_use, rownames(feat_mat))
#   missing <- setdiff(features_use, rownames(feat_mat))
#
#   if (length(present) < 1) {
#     stop("None of the requested features are present in the assay rownames().")
#   }
#   if (length(missing) > 0) {
#     warning("Dropping missing features (not found in rownames(assay)): ",
#             paste(head(missing, 20), collapse = ", "),
#             if (length(missing) > 20) paste0(" ... (+", length(missing) - 20, " more)") else "")
#   }
#
#   feat_mat <- feat_mat[present, , drop = FALSE]
#
#   feature_df <- as.data.frame(t(feat_mat))
#   feature_df$cell_id <- colnames(spe)
#
#   df <- feature_df
#   df[[label_col]] <- lab_vec
#
#   labelled_df <- df[!(is.na(df[[label_col]]) | df[[label_col]] == unknown_label), , drop = FALSE]
#
#   if (nrow(labelled_df) < 2) stop("Not enough labelled cells to train model.")
#
#   labelled_df[[label_col]] <- factor(labelled_df[[label_col]])
#   n_classes <- nlevels(labelled_df[[label_col]])
#   if (n_classes < 2) stop("At least two classes are required to train the model.")
#
#   if (!is.null(seed)) set.seed(seed)
#
#   p <- ncol(labelled_df) - 2  # exclude cell_id and label
#   if (p < 1) stop("Feature matrix has zero columns; cannot train model.")
#   mtry_val <- if (is.null(mtry)) max(1, floor(sqrt(p))) else mtry
#
#   # Cross-validation always (degenerate case cv_folds = 1 allowed)
#   cv_folds <- max(1, min(cv_folds, nrow(labelled_df)))
#   fold_assign <- integer(nrow(labelled_df))
#   for (cls in levels(labelled_df[[label_col]])) {
#     idx <- which(labelled_df[[label_col]] == cls)
#     fold_ids <- rep(seq_len(cv_folds), length.out = length(idx))
#     fold_assign[idx] <- sample(fold_ids)
#   }
#
#   formula <- stats::as.formula(paste(label_col, "~ ."))
#
#   f1_macro <- function(truth, pred) {
#     classes <- union(levels(truth), levels(pred))
#     f1s <- sapply(classes, function(cls) {
#       tp <- sum(truth == cls & pred == cls)
#       fp <- sum(truth != cls & pred == cls)
#       fn <- sum(truth == cls & pred != cls)
#       prec <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
#       rec <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
#       ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
#     })
#     mean(f1s)
#   }
#
#   overall_list <- vector("list", cv_folds)
#   class_list <- vector("list", cv_folds)
#   all_truth <- labelled_df[[label_col]]
#   all_pred <- factor(rep(NA_character_, nrow(labelled_df)), levels = levels(all_truth))
#
#   for (k in seq_len(cv_folds)) {
#     train_k <- labelled_df[fold_assign != k, , drop = FALSE]
#     test_k  <- labelled_df[fold_assign == k, , drop = FALSE]
#
#     if (nrow(train_k) < 2 || nrow(test_k) < 1) next
#
#     train_k$cell_id <- NULL
#     test_k$cell_id <- NULL
#
#     model_k <- ranger::ranger(
#       formula,
#       data = train_k,
#       num.trees = num.trees,
#       mtry = mtry_val,
#       num.threads = max(1, parallel::detectCores(logical = FALSE))
#     )
#
#     pred_k <- stats::predict(model_k, test_k)$predictions
#     pred_k <- factor(pred_k, levels = levels(all_truth))
#     truth_k <- test_k[[label_col]]
#
#     all_pred[fold_assign == k] <- pred_k
#
#     overall_list[[k]] <- data.frame(
#       fold = k,
#       accuracy = mean(pred_k == truth_k),
#       f1_macro = f1_macro(truth_k, pred_k)
#     )
#
#     classes <- union(levels(truth_k), levels(pred_k))
#     class_list[[k]] <- data.frame(
#       class = factor(classes, levels = classes),
#       tp = sapply(classes, function(cls) sum(truth_k == cls & pred_k == cls)),
#       fp = sapply(classes, function(cls) sum(truth_k != cls & pred_k == cls)),
#       fn = sapply(classes, function(cls) sum(truth_k == cls & pred_k != cls)),
#       support = sapply(classes, function(cls) sum(truth_k == cls)),
#       fold = k
#     ) |>
#       dplyr::mutate(
#         precision = ifelse(.data$tp + .data$fp == 0, 0, .data$tp / (.data$tp + .data$fp)),
#         recall    = ifelse(.data$tp + .data$fn == 0, 0, .data$tp / (.data$tp + .data$fn)),
#         f1        = ifelse(.data$precision + .data$recall == 0, 0,
#                            2 * .data$precision * .data$recall / (.data$precision + .data$recall))
#       )
#   }
#
#   overall_list <- Filter(Negate(is.null), overall_list)
#   class_list <- Filter(Negate(is.null), class_list)
#   cv_overall <- if (length(overall_list)) dplyr::bind_rows(overall_list) else NULL
#   cv_class_metrics <- if (length(class_list)) dplyr::bind_rows(class_list) else NULL
#
#   # Aggregate CV predictions for summary metrics
#   confusion_cv <- table(truth = all_truth, pred = all_pred)
#   acc_cv <- mean(all_truth == all_pred, na.rm = TRUE)
#   f1_cv <- if (any(is.na(all_pred))) NA_real_ else f1_macro(all_truth, all_pred)
#
#   metrics <- list(
#     confusion_matrix = confusion_cv,
#     accuracy = acc_cv,
#     f1_macro = f1_cv,
#     cv_overall = cv_overall
#   )
#
#   # Train final model on all labelled data
#   train_all <- labelled_df
#   train_all$cell_id <- NULL
#   final_model <- ranger::ranger(
#     formula,
#     data = train_all,
#     num.trees = num.trees,
#     mtry = mtry_val,
#     num.threads = max(1, parallel::detectCores(logical = FALSE))
#   )
#
#   list(
#     model = final_model,
#     metrics = metrics,
#     test_pred = all_pred,
#     test_truth = all_truth,
#     cv_overall = cv_overall,
#     cv_class_metrics = cv_class_metrics,
#     features_used = present   # NEW: return what was actually used
#   )
# }
# train_custom_randomforest <- function(spe,
#                                       label_col = "cutoff_label",
#                                       assay_name = "norm",
#                                       unknown_label = "Unknown",
#                                       num.trees = 200,
#                                       mtry = NULL,
#                                       seed = NULL,
#                                       features = "all",
#                                       cv_folds = 5,
#                                       repeats = 10,
#                                       agreement_thresh = 0.8) {
#
#   .assert_spe(spe)
#
#   # 1. Feature Selection Logic
#   feat_mat_all <- SummarizedExperiment::assay(spe, assay_name)
#
#   if (is.character(features)) {
#     if (length(features) == 1 && features == "all") {
#       features_use <- rownames(feat_mat_all)
#     } else {
#       features_use <- intersect(features, rownames(feat_mat_all))
#     }
#   } else {
#     stop("features must be 'all' or a character vector of markers.")
#   }
#
#   if (length(features_use) < 1) stop("No valid features selected.")
#
#   # 2. Prepare Data
#   lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
#   feature_df <- as.data.frame(t(feat_mat_all[features_use, , drop = FALSE]))
#   feature_df$original_label <- lab_vec
#
#   # Identify index of core cells (those not labeled Unknown)
#   core_idx <- which(feature_df$original_label != unknown_label & !is.na(feature_df$original_label))
#   if (length(core_idx) < 10) stop("Not enough core cells for CV cleaning.")
#
#   core_df <- feature_df[core_idx, ]
#   core_df$original_label <- factor(core_df$original_label)
#
#   # 3. Repeated k-Fold CV for Label Cleaning
#   if (!is.null(seed)) set.seed(seed)
#   match_counts <- setNames(numeric(nrow(core_df)), rownames(core_df))
#
#   message(sprintf("Starting label cleaning: %d repeats of %d-fold CV...", repeats, cv_folds))
#
#   for (r in seq_len(repeats)) {
#     message(sprintf("Repeat: %d",r))
#     fold_assign <- integer(nrow(core_df))
#     for (cls in levels(core_df$original_label)) {
#       cls_idx <- which(core_df$original_label == cls)
#       fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
#     }
#
#     for (k in seq_len(cv_folds)) {
#       train_idx <- which(fold_assign != k)
#       test_idx  <- which(fold_assign == k)
#
#       temp_model <- ranger::ranger(
#         original_label ~ .,
#         data = core_df[train_idx, ],
#         num.trees = 100,
#         num.threads = max(1, parallel::detectCores(logical = FALSE))
#       )
#
#       preds <- stats::predict(temp_model, core_df[test_idx, ])$predictions
#       is_match <- (preds == core_df$original_label[test_idx])
#       match_counts[test_idx] <- match_counts[test_idx] + as.numeric(is_match)
#     }
#   }
#
#   # 4. Filter Core Cells based on Consensus
#   agreement_rate <- match_counts / repeats
#   # Logic: Keep if agreement >= threshold
#   valid_core_names <- names(agreement_rate)[agreement_rate >= agreement_thresh]
#   inconsistent_names <- names(agreement_rate)[agreement_rate < agreement_thresh]
#
#   # Update labels in the SPE object
#   cleaned_labels <- lab_vec
#   # Match by cell names (colnames of SPE)
#   cleaned_labels[colnames(spe) %in% inconsistent_names] <- unknown_label
#   SummarizedExperiment::colData(spe)$cleaned_core_label <- cleaned_labels
#
#   message(sprintf("Cleaning complete. Removed %d inconsistent core cells.",
#                   length(inconsistent_names)))
#
#   # 5. Final Model Training on Cleaned Set
#   final_train_df <- core_df[valid_core_names, ]
#
#   # Handle mtry default if NULL
#   mtry_val <- if (is.null(mtry)) floor(sqrt(length(features_use))) else mtry
#
#   final_model <- ranger::ranger(
#     original_label ~ .,
#     data = final_train_df,
#     num.trees = num.trees,
#     mtry = mtry_val,probability = TRUE,
#     num.threads = max(1, parallel::detectCores(logical = FALSE))
#   )
#
#   # 6. Reporting and Return
#   return(list(
#     spe = spe,
#     model = final_model,
#     features_used = features_use,
#     agreement_rates = agreement_rate,
#     cleaned_core_names = valid_core_names
#   ))
# }
#
# train_custom_randomforest <- function(spe,
#                                       label_col = "cutoff_label",
#                                       assay_name = "norm",
#                                       unknown_label = "Unknown",
#                                       num.trees = 200,
#                                       mtry = NULL,
#                                       seed = NULL,
#                                       features = "all",
#                                       cv_folds = 5,
#                                       repeats = 10,
#                                       agreement_thresh = 0.8,
#                                       # parallel = TRUE,
#                                       num_threads=2) {
#
#   if (!is.null(seed)) set.seed(seed)
#   .assert_spe(spe)
#
#   # num_threads <- if (parallel) max(1, parallel::detectCores(logical = FALSE)) else 1
#
#   # 1. Feature Prep
#   feat_mat_all <- SummarizedExperiment::assay(spe, assay_name)
#   features_use <- if (length(features) == 1 && features == "all") rownames(feat_mat_all) else intersect(features, rownames(feat_mat_all))
#
#   lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
#   feature_df <- as.data.frame(t(feat_mat_all[features_use, , drop = FALSE]))
#
#   # Filter only labelled cells for training/cleaning
#   core_idx <- which(lab_vec != unknown_label & !is.na(lab_vec))
#   core_df <- feature_df[core_idx, ]
#   core_df$original_label <- factor(lab_vec[core_idx])
#
#   # --- HELPER: F1 MACRO ---
#   f1_macro <- function(truth, pred) {
#     classes <- levels(truth)
#     f1s <- sapply(classes, function(cls) {
#       tp <- sum(truth == cls & pred == cls)
#       fp <- sum(truth != cls & pred == cls)
#       fn <- sum(truth == cls & pred != cls)
#       prec <- if(tp + fp == 0) 0 else tp / (tp + fp)
#       rec  <- if(tp + fn == 0) 0 else tp / (tp + fn)
#       if(prec + rec == 0) 0 else 2 * (prec * rec) / (prec + rec)
#     })
#     return(f1s) # Returns vector for all classes
#   }
#
#   # --- STAGE 1: CONSENSUS CLEANING LOOP ---
#   message(sprintf("Stage 1: Cleaning labels via %d repeats...", repeats))
#   match_counts <- setNames(numeric(nrow(core_df)), rownames(core_df))
#
#   for (r in seq_len(repeats)) {
#     fold_assign <- integer(nrow(core_df))
#     for (cls in levels(core_df$original_label)) {
#       cls_idx <- which(core_df$original_label == cls)
#       fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
#     }
#     for (k in seq_len(cv_folds)) {
#       train_idx <- which(fold_assign != k); test_idx <- which(fold_assign == k)
#       tmp <- ranger::ranger(original_label ~ ., data = core_df[train_idx, ], num.trees = 100, num.threads = num_threads)
#       preds <- stats::predict(tmp, core_df[test_idx, ])$predictions
#       match_counts[test_idx] <- match_counts[test_idx] + as.numeric(preds == core_df$original_label[test_idx])
#     }
#   }
#
#   agreement_rate <- match_counts / repeats
#   valid_names <- names(agreement_rate)[agreement_rate >= agreement_thresh]
#   cleaned_df <- core_df[valid_names, ]
#
#   # --- STAGE 2: EVALUATION LOOP (On Cleaned Data Only) ---
#   message("Stage 2: Evaluating final model performance on cleaned cells...")
#   eval_preds <- factor(rep(NA_character_, nrow(cleaned_df)), levels = levels(cleaned_df$original_label))
#   eval_fold_assign <- integer(nrow(cleaned_df))
#   for (cls in levels(cleaned_df$original_label)) {
#     cls_idx <- which(cleaned_df$original_label == cls)
#     eval_fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
#   }
#
#   overall_list <- list()
#   class_metric_list <- list()
#
#
#
#   for (k in seq_len(cv_folds)) {
#     train_idx <- which(eval_fold_assign != k); test_idx <- which(eval_fold_assign == k)
#     eval_mod <- ranger::ranger(original_label ~ ., data = cleaned_df[train_idx, ], num.trees = num.trees, num.threads = num_threads)
#     pk <- stats::predict(eval_mod, cleaned_df[test_idx, ])$predictions
#     pk <- factor(pk, levels = levels(cleaned_df$original_label))
#     tk <- cleaned_df$original_label[test_idx]
#
#     eval_preds[test_idx] <- pk
#
#     # Calculate Overall Fold Metrics
#     f1_vec <- f1_macro(tk, pk)
#     overall_list[[k]] <- data.frame(fold = k, accuracy = mean(pk == tk), f1_macro = mean(f1_vec))
#
#     # Calculate Per-Class Metrics for this fold
#     classes <- levels(cleaned_df$original_label)
#     class_metric_list[[k]] <- data.frame(
#       class = classes,
#       fold = k,
#       precision = sapply(classes, function(c) {
#         tp <- sum(tk==c & pk==c); fp <- sum(tk!=c & pk==c)
#         if(tp+fp==0) 0 else tp/(tp+fp)
#       }),
#       recall = sapply(classes, function(c) {
#         tp <- sum(tk==c & pk==c); fn <- sum(tk==c & pk!=c)
#         if(tp+fn==0) 0 else tp/(tp+fn)
#       }),
#       f1 = f1_vec
#     )
#   }
#
#   # --- FINAL MODEL TRAINING ---
#   mtry_val <- if (is.null(mtry)) floor(sqrt(length(features_use))) else mtry
#   final_model <- ranger::ranger(original_label ~ ., data = cleaned_df,
#                                 num.trees = num.trees, mtry = mtry_val,
#                                 probability = TRUE, num.threads = num_threads)
#
#   # Update SPE colData
#   new_labels <- lab_vec
#   new_labels[colnames(spe) %in% names(agreement_rate[agreement_rate < agreement_thresh])] <- unknown_label
#   SummarizedExperiment::colData(spe)$cleaned_core_label <- new_labels
#
#   # --- RETURN LIST (Matches your old output structure) ---
#   return(list(
#     spe = spe,
#     model = final_model,
#     agreement_rates = agreement_rate,
#     features_used = features_use,
#     metrics = list(
#       confusion_matrix = table(truth = cleaned_df$original_label, pred = eval_preds),
#       accuracy = mean(cleaned_df$original_label == eval_preds),
#       f1_macro = mean(f1_macro(cleaned_df$original_label, eval_preds)),
#       cv_overall = dplyr::bind_rows(overall_list),
#       cv_class_metrics = dplyr::bind_rows(class_metric_list)
#     ),
#     test_pred = eval_preds,
#     test_truth = cleaned_df$original_label
#   ))
# }



#' Predict and fill unknown labels using a trained random forest model
#'
#' Use a trained \code{ranger} classification model to predict cell-type labels
#' for all cells in a \code{SpatialExperiment} or \code{SingleCellExperiment}
#' object. Predicted labels are used to fill entries in \code{label_col} that
#' are either \code{NA} or equal to \code{unknown_label}.
#'
#' If the model was trained with \code{probability = TRUE}, class probabilities
#' are used to compute a confidence score for each prediction. Predictions with
#' maximum class probability below \code{threshold} are reassigned to
#' \code{unassigned_label}. The maximum probability is stored for quality control.
#'
#' @param spe A \code{SpatialExperiment} or \code{SingleCellExperiment} object.
#' @param model A fitted \code{ranger} model returned by \code{train_custom_rf()}.
#'   The model should be trained with \code{probability = TRUE} to enable
#'   threshold-based filtering.
#' @param assay_name Character scalar specifying which assay to use as feature
#'   input (default \code{"norm"}).
#' @param label_col Character scalar naming the existing label column in
#'   \code{colData(spe)} (default \code{"custom_label"}).
#' @param out_col Character scalar naming the output column that will contain
#'   the filled labels (default \code{"soft_tree_label_filled"}).
#' @param pred_col Character scalar naming the column used to store raw model
#'   predictions (default \code{"rf_pred"}).
#' @param unknown_label Character value treated as missing and eligible for
#'   replacement (default \code{"Unknown"}).
#' @param threshold Numeric value in [0,1] specifying the minimum required
#'   class probability for a prediction to be accepted (default \code{0.5}).
#'   Ignored if the model was not trained with \code{probability = TRUE}.
#' @param unassigned_label Character label assigned to cells whose maximum
#'   class probability is below \code{threshold} (default \code{"Unassigned"}).
#'
#' @return The input \code{spe} object with updated \code{colData} columns:
#'   \describe{
#'     \item{\code{pred_col}}{Raw predicted labels from the random forest.}
#'     \item{\code{out_col}}{Filled labels, replacing only \code{NA} or
#'       \code{unknown_label} entries in \code{label_col}.}
#'     \item{\code{rf_confidence}}{Maximum predicted class probability per cell.}
#'   }
#'
#' @details
#' Feature columns are automatically aligned to match the variables used during
#' model training. Only cells originally labeled as \code{NA} or
#' \code{unknown_label} are replaced in the output column.
#'
#' @export
# predict_unknown_with_rf <- function(spe,
#                                     model,
#                                     assay_name = "norm",
#                                     label_col = "custom_label",
#                                     out_col = "soft_tree_label_filled",
#                                     pred_col = "rf_pred",
#                                     unknown_label = "Unknown") {
#   .assert_spe(spe)
#   if (!inherits(model, "ranger")) stop("model must be a ranger object.")
#   if (!label_col %in% colnames(SummarizedExperiment::colData(spe))) {
#     stop("label_col not found in colData(spe).")
#   }
#   if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
#     stop(paste0("Assay '", assay_name, "' not found in spe."))
#   }
#
#   feat_mat <- SummarizedExperiment::assay(spe, assay_name)
#   feature_df <- as.data.frame(t(feat_mat))
#
#   # Align feature columns with the model's expected variables
#   expected <- model$forest$independent.variable.names
#   missing_cols <- setdiff(expected, names(feature_df))
#   if (length(missing_cols) > 0) {
#     stop("Missing feature columns required by model: ", paste(missing_cols, collapse = ", "))
#   }
#   feature_df <- feature_df[, expected, drop = FALSE]
#
#   preds <- stats::predict(model, feature_df)$predictions
#
#   labels <- SummarizedExperiment::colData(spe)[[label_col]]
#   filled <- labels
#   replace_idx <- is.na(labels) | labels == unknown_label
#   filled[replace_idx] <- as.character(preds)[replace_idx]
#
#   SummarizedExperiment::colData(spe)[[pred_col]] <- preds
#   SummarizedExperiment::colData(spe)[[out_col]] <- filled
#   spe
# }


###############BIC



#' Convert a tree to a tibble representation
#'
#' @param node A gating tree (node or leaf).
#' @param parent Parent node id (used internally for recursion).
#' @param id Node id for the current subtree (used internally for recursion).
#'
#' @return A tibble with columns `id`, `parent`, `type`, `label`.
#' @export

tree_to_df <- function(node, parent = NA_character_, id = "root") {
  if (node$type == "leaf") {
    return(tibble::tibble(
      id = id,
      parent = parent,
      type = "leaf",
      label = paste0("leaf\nn=", length(node$cells))
    ))
  }

  this <- tibble::tibble(
    id = id,
    parent = parent,
    type = "node",
    label = sprintf("%s > %.2f\nsep=%.2f  w=%.2f",
                    node$marker,
                    node$cutoff,
                    node$sep_score %||% NA_real_,
                    node$w %||% NA_real_)
  )

  left <- tree_to_df(node$left, id, paste0(id, "_L"))
  right <- tree_to_df(node$right, id, paste0(id, "_R"))

  dplyr::bind_rows(this, left, right)
}



#' Print a cell-type gating tree
#'
#' @param node A gating tree (node or leaf) produced by [build_fullcoverage_tree()].
#' @param indent String used internally for indentation during recursive printing.
#'
#' @return Invisibly returns `NULL`.
#' @export
print_celltype_tree <- function(node, indent = "") {
  if (node$type == "leaf") {
    cat(indent, "[leaf] depth =", node$depth,
        " n_cells =", length(node$cells), "\n", sep = "")
    return(invisible(NULL))
  }

  cat(indent,
      sprintf("[depth %d] %s > %.3f  (sep = %.2f, w = %.2f, n = %d)\n",
              node$depth,
              node$marker,
              node$cutoff,
              node$sep_score %||% NA_real_,
              node$w %||% NA_real_,
              length(node$cells)),
      sep = "")

  cat(indent, " |- low\n", sep = "")
  print_celltype_tree(node$left, paste0(indent, " |  "))

  cat(indent, " `- high\n", sep = "")
  print_celltype_tree(node$right, paste0(indent, "    "))
}



#' Fit a 2-component Gaussian mixture model
#'
#' @param x numeric vector
#' @param cutoff_method Character. One of "mean" (midpoint of means) or
#'   "equal_posteriors" (threshold where component posteriors are equal).
#' @param gmm_model_names Optional character vector of model names to pass to
#'   [mclust::Mclust()] (e.g., "V" to forbid equal-variance in 1D).
#' @return list with GMM parameters or NULL
#' @export
fit_gmm_2 <- function(x,
                      cutoff_method = c("mean", "equal_posteriors"),
                      gmm_model_names = NULL) {
  cutoff_method <- match.arg(cutoff_method)
  x <- x[is.finite(x)]
  if (length(x) < 50L || length(unique(x)) < 3L) return(NULL)

  gmm2 <- tryCatch(
    mclust::Mclust(x, G = 2, modelNames = gmm_model_names, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(gmm2) || length(gmm2$parameters$mean) != 2) return(NULL)

  ## NEW: fit the unimodal (G=1) competitor and compare via BIC.
  ## mclust's BIC convention is HIGHER = better fit.
  gmm1 <- tryCatch(
    mclust::Mclust(x, G = 1, modelNames = gmm_model_names, verbose = FALSE),
    error = function(e) NULL
  )

  bic2 <- gmm2$bic
  bic1 <- if (!is.null(gmm1)) gmm1$bic else NA_real_

  ## NEW: BIC-to-posterior-probability approximation (Kass & Raftery 1995).
  ## w -> 1 means strong evidence for two components; w -> 0.5 means the
  ## two models are indistinguishable; w below 0.5 would mean G=1 is
  ## actually favored (shouldn't normally happen since we picked gmm2,
  ## but w still reports the true relative evidence either way).
  if (is.finite(bic1) && is.finite(bic2)) {
    m <- max(bic1, bic2)  # subtract max for numerical stability before exponentiating
    w <- exp((bic2 - m) / 2) / (exp((bic2 - m) / 2) + exp((bic1 - m) / 2))
  } else {
    ## G=1 fit failed to converge -- no competing model to compare against,
    ## so there's no evidence AGAINST bimodality. Default to full trust.
    w <- 1
  }

  mu <- as.numeric(gmm2$parameters$mean)

  sig2 <- gmm2$parameters$variance$sigmasq
  if (length(sig2) == 1) {
    sig2 <- rep(sig2, 2)
  }
  sig2 <- as.numeric(sig2)
  pi <- as.numeric(gmm2$parameters$pro)

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
       cutoff = cutoff, sep_score = sep,
       w = w, bic1 = bic1, bic2 = bic2)   ## NEW fields
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
#' @export

build_fullcoverage_tree <- function(expr_mat,
                                    markers_pos,
                                    markers_neg = character(0),
                                    cell_idx = seq_len(ncol(expr_mat)),
                                    depth = 0,
                                    max_depth = 4,
                                    min_cells = 200,
                                    min_score = 0.5,
                                    w_threshold = 0.5,   ## NEW argument
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
    fit <- fit_gmm_2(x, cutoff_method = cutoff_method, gmm_model_names = gmm_model_names)
    if (is.null(fit)) {
      return(data.frame(marker = m, sep = NA_real_, acc = NA_real_,
                        cutoff = NA_real_, scale = NA_real_, w = NA_real_))  ## NEW: w column
    }

    cutoff <- fit$cutoff
    scale <- sqrt(fit$s1^2 + fit$s2^2)
    if (!is.finite(scale) || scale == 0) scale <- .safe_mad(x)

    pred <- ifelse(x > cutoff, 1, 0)
    acc <- max(mean(pred == target, na.rm = TRUE),
               mean((1 - pred) == target, na.rm = TRUE))

    data.frame(marker = m, sep = fit$sep_score, acc = acc,
               cutoff = cutoff, scale = scale, w = fit$w)  ## NEW: carry w through
  })
  print(marker_stats)
  marker_stats <- do.call(rbind, marker_stats)
  marker_stats$sep <- as.numeric(marker_stats$sep)
  marker_stats$acc <- as.numeric(marker_stats$acc)
  marker_stats$w   <- as.numeric(marker_stats$w)          ## NEW
  marker_stats$combo <- marker_stats$sep + 2 * (marker_stats$acc - 0.5)

  marker_stats <- marker_stats[order(-marker_stats$combo), , drop = FALSE]

  ## NEW: eligibility scan. Ranking is unchanged (still combo-based); the
  ## difference is we no longer only look at row 1 -- we walk down the
  ## combo-sorted list until we find the first marker that ALSO clears the
  ## bimodality-evidence bar. A marginal top-combo-but-low-w marker no
  ## longer blocks a solid, genuinely bimodal runner-up.
  eligible <- which(
    is.finite(marker_stats$sep) & marker_stats$sep >= min_score &
      is.finite(marker_stats$w)   & marker_stats$w   >= w_threshold
  )

  if (length(eligible) == 0) {
    return(list(type = "leaf", depth = depth, cells = cell_idx))
  }

  best <- marker_stats[eligible[1], ]
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
    w = as.numeric(best$w),                                ## NEW: stored on the node
    cells = cell_idx,
    left = build_fullcoverage_tree(expr_mat, remaining_pos, markers_neg, left,
                                   depth + 1, max_depth, min_cells, min_score, w_threshold,
                                   cutoff_method = cutoff_method,
                                   gmm_model_names = gmm_model_names),
    right = build_fullcoverage_tree(expr_mat, remaining_pos, markers_neg, right,
                                    depth + 1, max_depth, min_cells, min_score, w_threshold,
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
#' @export
collect_path_scores <- function(tree, expr_mat, cell_i) {
  scores <- numeric(0)
  node <- tree
  while (!is.null(node) && node$type == "node") {
    m <- node$marker
    x <- as.numeric(expr_mat[m, cell_i])
    raw <- score_marker_logistic(x, node$cutoff, node$scale)

    ## NEW: shrink the raw score toward the uninformative value (0.5) in
    ## proportion to how much bimodal evidence this node's marker actually
    ## had. w = 1 -> full trust, raw score unchanged. w close to
    ## w_threshold -> score pulled most of the way to 0.5, barely moving
    ## the combined probability either direction.
    w <- node$w %||% 1  ## backward-compat: trees built before this change
    ## have no $w field, default to full trust
    s <- w * raw + (1 - w) * 0.5

    scores <- c(scores, s)
    if (is.finite(x) && x > node$cutoff) node <- node$right else node <- node$left
  }
  scores
}
