#' Flag confident cells from a probability matrix
#'
#' Compute per-cell confidence using a user-supplied decision function or a
#' robust Tukey-style rule: a cell is considered confident if its maximum
#' lineage probability exceeds `max(base_thresh, quantile + iqr_mult * IQR)` for
#' that cell.
#'
#' @param prob_mat Numeric matrix of probabilities (cells x lineages).
#' @param base_thresh Minimum allowable threshold (default 0.4) used only for the
#'   default rule.
#' @param quantile_cut Quantile used to set the adaptive threshold (default 0.75)
#'   for the default rule.
#' @param iqr_mult Multiplier applied to the IQR (default 1.5) for the default
#'   rule.
#' @param flag_fn Optional function taking a numeric vector of probabilities and
#'   returning a single logical value. Defaults to the Tukey-style rule above.
#'
#' @return Logical vector length `nrow(prob_mat)` indicating confident cells.
#' @export
compute_custom_labels <- function(prob_mat,
                                   base_thresh = 0.4,
                                   quantile_cut = 0.75,
                                   iqr_mult = 1.5,
                                   flag_fn = NULL) {
  if (!is.matrix(prob_mat)) stop("prob_mat must be a matrix.")
  if (nrow(prob_mat) == 0 || ncol(prob_mat) == 0) {
    stop("prob_mat must have at least one row and one column.")
  }

  default_flag <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0) return(FALSE)
    q <- stats::quantile(x, quantile_cut, na.rm = TRUE, names = FALSE)
    iqr <- stats::IQR(x, na.rm = TRUE)
    thr <- max(base_thresh, q + iqr_mult * iqr)
    max(x, na.rm = TRUE) > thr
  }

  fn <- flag_fn %||% default_flag

  flags <- apply(prob_mat, 1, function(x) {
    out <- fn(x)
    if (length(out) != 1 || !is.logical(out)) {
      stop("flag_fn must return a single logical value per cell.")
    }
    isTRUE(out)
  })

  as.logical(flags)
}


#' Assign labels with uncertainty handling
#'
#' Use a probability matrix and confidence flags to assign lineage labels;
#' cells that are not confident receive `unknown_label`.
#'
#' @param prob_mat Numeric matrix of probabilities (cells x lineages).
#' @param confidence_flags Logical vector length `nrow(prob_mat)`.
#' @param unknown_label Character label for uncertain cells (default "Unknown").
#'
#' @return Character vector of labels length `nrow(prob_mat)`.
#' @export
assign_confident_labels <- function(prob_mat,
                                    confidence_flags,
                                    unknown_label = "Unknown") {
  if (!is.matrix(prob_mat)) stop("prob_mat must be a matrix.")
  if (nrow(prob_mat) == 0 || ncol(prob_mat) == 0) {
    stop("prob_mat must have at least one row and one column.")
  }
  if (length(confidence_flags) != nrow(prob_mat)) {
    stop("confidence_flags must have length equal to nrow(prob_mat).")
  }

  max_probs <- apply(prob_mat, 1, function(x) {
    x <- as.numeric(x)
    if (all(!is.finite(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  })

  lab <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
  lab[!is.finite(max_probs)] <- unknown_label
  lab[!confidence_flags] <- unknown_label
  lab
}


#' Add custom labels to a SpatialExperiment
#'
#' Applies the confidence-based labelling scheme to a probability matrix and
#' stores the resulting labels in `colData(spe)`.
#'
#' @param spe A `SpatialExperiment` or `SingleCellExperiment`.
#' @param prob_mat Numeric probability matrix (cells x lineages). Rows must
#'   correspond to `colnames(spe)`; if rownames are present they are checked.
#' @param colname Name of the output column in `colData` (default "custom_label").
#' @param base_thresh Minimum allowable threshold (default 0.4).
#' @param quantile_cut Quantile used to set the adaptive threshold (default 0.75).
#' @param iqr_mult Multiplier applied to the IQR (default 1.5).
#' @param unknown_label Label assigned to uncertain cells (default "Unknown").
#' @param flag_fn Optional function taking a probability vector and returning a
#'   logical indicating confidence; defaults to the internal Tukey-style rule.
#'
#' @return The input `spe` with a new `colData` column containing custom labels.
#' @export
custom_labels <- function(spe,
                          prob_mat,
                          colname = "custom_label",
                          base_thresh = 0.4,
                          quantile_cut = 0.75,
                          iqr_mult = 1.5,
                          unknown_label = "Unknown",
                          flag_fn = NULL) {
  .assert_spe(spe)
  if (!is.matrix(prob_mat)) stop("prob_mat must be a matrix.")
  if (nrow(prob_mat) != ncol(spe)) {
    stop("nrow(prob_mat) must match ncol(spe).")
  }
  if (!is.null(rownames(prob_mat))) {
    if (!identical(rownames(prob_mat), colnames(spe))) {
      stop("rownames(prob_mat) must match colnames(spe) when provided.")
    }
  }

  flags <- compute_custom_labels(prob_mat,
                                 base_thresh = base_thresh,
                                 quantile_cut = quantile_cut,
                                 iqr_mult = iqr_mult,
                                 flag_fn = flag_fn)
  labels <- assign_confident_labels(prob_mat, flags, unknown_label = unknown_label)

  SummarizedExperiment::colData(spe)[[colname]] <- labels
  spe
}


#' Train a random forest on labelled cells
#'
#' Trains a `ranger` random forest using assay features and an existing label
#' column in `colData`. Cells with missing or `unknown_label` are excluded from
#' training and evaluation.
#'
#' @param spe A `SpatialExperiment` or `SingleCellExperiment`.
#' @param label_col Name of the label column in `colData` (default "custom_label").
#' @param assay_name Assay to use as features (default "norm").
#' @param unknown_label Label value to treat as unlabeled (default "Unknown").
#' @param train_frac Fraction of labelled cells used for training (currently
#'   ignored; cross-validation uses all labelled cells).
#' @param num.trees Number of trees for `ranger` (default 200).
#' @param mtry Optional `mtry`; if `NULL`, uses floor(sqrt(p)).
#' @param seed Optional seed for reproducibility.
#' @param cv_folds Number of folds for cross-validation; if >1, CV metrics are
#'   computed in addition to the hold-out split (default 5).
#'
#' @return List with `model`, `metrics` (CV confusion_matrix, accuracy,
#'   f1_macro, cv_overall), `test_pred`/`test_truth` (pooled CV predictions and
#'   truths), and `cv_overall`/`cv_class_metrics` per fold.
#' @export
train_custom_rf <- function(spe,
                            label_col = "custom_label",
                            assay_name = "norm",
                            unknown_label = "Unknown",
                            train_frac = 0.8,
                            num.trees = 200,
                            mtry = NULL,
                            seed = NULL,
                            cv_folds = 5) {
  .assert_spe(spe)
  if (!label_col %in% colnames(SummarizedExperiment::colData(spe))) {
    stop("label_col not found in colData(spe).")
  }
  if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
    stop(paste0("Assay '", assay_name, "' not found in spe."))
  }

  lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
  feat_mat <- SummarizedExperiment::assay(spe, assay_name)

  feature_df <- as.data.frame(t(feat_mat))
  feature_df$cell_id <- colnames(spe)

  df <- feature_df
  df[[label_col]] <- lab_vec

  labelled_df <- df[!(is.na(df[[label_col]]) | df[[label_col]] == unknown_label), , drop = FALSE]

  if (nrow(labelled_df) < 2) stop("Not enough labelled cells to train model.")

  labelled_df[[label_col]] <- factor(labelled_df[[label_col]])
  n_classes <- nlevels(labelled_df[[label_col]])
  if (n_classes < 2) stop("At least two classes are required to train the model.")

  if (!is.null(seed)) set.seed(seed)

  p <- ncol(labelled_df) - 2  # exclude cell_id and label
  if (p < 1) stop("Feature matrix has zero columns; cannot train model.")
  mtry_val <- if (is.null(mtry)) max(1, floor(sqrt(p))) else mtry

  # Cross-validation always (degenerate case cv_folds = 1 allowed)
  cv_folds <- max(1, min(cv_folds, nrow(labelled_df)))
  fold_assign <- integer(nrow(labelled_df))
  for (cls in levels(labelled_df[[label_col]])) {
    idx <- which(labelled_df[[label_col]] == cls)
    fold_ids <- rep(seq_len(cv_folds), length.out = length(idx))
    fold_assign[idx] <- sample(fold_ids)
  }

  formula <- stats::as.formula(paste(label_col, "~ ."))

  f1_macro <- function(truth, pred) {
    classes <- union(levels(truth), levels(pred))
    f1s <- sapply(classes, function(cls) {
      tp <- sum(truth == cls & pred == cls)
      fp <- sum(truth != cls & pred == cls)
      fn <- sum(truth == cls & pred != cls)
      prec <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
      rec <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
      ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
    })
    mean(f1s)
  }

  overall_list <- vector("list", cv_folds)
  class_list <- vector("list", cv_folds)
  all_truth <- labelled_df[[label_col]]
  all_pred <- factor(rep(NA_character_, nrow(labelled_df)), levels = levels(all_truth))

  for (k in seq_len(cv_folds)) {
    train_k <- labelled_df[fold_assign != k, , drop = FALSE]
    test_k  <- labelled_df[fold_assign == k, , drop = FALSE]

    if (nrow(train_k) < 2 || nrow(test_k) < 1) next

    train_k$cell_id <- NULL
    test_k$cell_id <- NULL

    model_k <- ranger::ranger(
      formula,
      data = train_k,
      num.trees = num.trees,
      mtry = mtry_val,
      num.threads = max(1, parallel::detectCores(logical = FALSE))
    )

    pred_k <- predict(model_k, test_k)$predictions
    pred_k <- factor(pred_k, levels = levels(all_truth))
    truth_k <- test_k[[label_col]]

    all_pred[fold_assign == k] <- pred_k

    overall_list[[k]] <- data.frame(
      fold = k,
      accuracy = mean(pred_k == truth_k),
      f1_macro = f1_macro(truth_k, pred_k)
    )

    classes <- union(levels(truth_k), levels(pred_k))
    class_list[[k]] <- data.frame(
      class = factor(classes, levels = classes),
      tp = sapply(classes, function(cls) sum(truth_k == cls & pred_k == cls)),
      fp = sapply(classes, function(cls) sum(truth_k != cls & pred_k == cls)),
      fn = sapply(classes, function(cls) sum(truth_k == cls & pred_k != cls)),
      support = sapply(classes, function(cls) sum(truth_k == cls)),
      fold = k
    ) |>
      dplyr::mutate(
        precision = ifelse(tp + fp == 0, 0, tp / (tp + fp)),
        recall = ifelse(tp + fn == 0, 0, tp / (tp + fn)),
        f1 = ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
      )
  }

  overall_list <- Filter(Negate(is.null), overall_list)
  class_list <- Filter(Negate(is.null), class_list)
  cv_overall <- if (length(overall_list)) dplyr::bind_rows(overall_list) else NULL
  cv_class_metrics <- if (length(class_list)) dplyr::bind_rows(class_list) else NULL

  # Aggregate CV predictions for summary metrics
  confusion_cv <- table(truth = all_truth, pred = all_pred)
  acc_cv <- mean(all_truth == all_pred, na.rm = TRUE)
  f1_cv <- if (any(is.na(all_pred))) NA_real_ else f1_macro(all_truth, all_pred)

  metrics <- list(
    confusion_matrix = confusion_cv,
    accuracy = acc_cv,
    f1_macro = f1_cv,
    cv_overall = cv_overall
  )

  # Train final model on all labelled data
  train_all <- labelled_df
  train_all$cell_id <- NULL
  final_model <- ranger::ranger(
    formula,
    data = train_all,
    num.trees = num.trees,
    mtry = mtry_val,
    num.threads = max(1, parallel::detectCores(logical = FALSE))
  )

  list(
    model = final_model,
    metrics = metrics,
    test_pred = all_pred,
    test_truth = all_truth,
    cv_overall = cv_overall,
    cv_class_metrics = cv_class_metrics
  )
}


#' Derive per-class metrics from a fitted model
#'
#' Prefers cross-validation class metrics when available; otherwise derives
#' metrics from the provided truth and prediction vectors.
#'
#' @param fit List returned by [train_custom_rf()].
#' @param test_truth Optional truth vector (used only if CV metrics absent).
#' @param test_pred Optional prediction vector (used only if CV metrics absent).
#'
#' @return Data frame with columns `class`, `tp`, `fp`, `fn`, `support`,
#'   `precision`, `recall`, `f1`, and optionally `fold` when sourced from CV.
#' @export
class_metrics_from_fit <- function(fit, test_truth = NULL, test_pred = NULL) {
  if (!is.null(fit$cv_class_metrics)) {
    return(fit$cv_class_metrics)
  }

  truth <- test_truth %||% fit$test_truth
  pred <- test_pred %||% fit$test_pred

  if (is.null(truth) || is.null(pred)) {
    stop("No class metrics available: provide test_truth/test_pred or ensure fit contains them.")
  }

  classes <- sort(unique(c(as.character(truth), as.character(pred))))
  out <- data.frame(
    class = factor(classes, levels = classes),
    tp = sapply(classes, function(cls) sum(truth == cls & pred == cls)),
    fp = sapply(classes, function(cls) sum(truth != cls & pred == cls)),
    fn = sapply(classes, function(cls) sum(truth == cls & pred != cls)),
    support = sapply(classes, function(cls) sum(truth == cls))
  ) |>
    dplyr::mutate(
      precision = ifelse(tp + fp == 0, 0, tp / (tp + fp)),
      recall = ifelse(tp + fn == 0, 0, tp / (tp + fn)),
      f1 = ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    )

  out
}


#' Confusion matrix between two label columns
#'
#' Builds a confusion matrix comparing two label columns in `colData(spe)`, with
#' aligned levels so missing categories are represented.
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment`.
#' @param label_col1 First label column name.
#' @param label_col2 Second label column name.
#' @param drop_na Logical; if TRUE, rows with NA in either label are dropped. Default TRUE.
#'
#' @return A table with rows = `label_col1` levels and columns = `label_col2` levels.
#' @export
label_confusion_matrix <- function(spe,
                                   label_col1,
                                   label_col2,
                                   drop_na = TRUE) {
  .assert_spe(spe)
  cd <- SummarizedExperiment::colData(spe)

  if (!label_col1 %in% colnames(cd)) {
    stop("label_col1 '", label_col1, "' not found in colData(spe).")
  }
  if (!label_col2 %in% colnames(cd)) {
    stop("label_col2 '", label_col2, "' not found in colData(spe).")
  }

  v1 <- cd[[label_col1]]
  v2 <- cd[[label_col2]]

  if (drop_na) {
    keep <- !(is.na(v1) | is.na(v2))
    v1 <- v1[keep]
    v2 <- v2[keep]
  }

  levels_all <- sort(unique(c(as.character(v1), as.character(v2))))
  v1 <- factor(v1, levels = levels_all)
  v2 <- factor(v2, levels = levels_all)

  table(v1, v2)
}


#' Agreement rates between two label columns
#'
#' Computes per-label agreement between two label columns in `colData(spe)`.
#' For each label (union of both columns), agreement is defined as the number of
#' cells where both labels match that value divided by the number of cells where
#' either column uses that value (Jaccard-style overlap).
#'
#' @param spe A `SpatialExperiment`/`SingleCellExperiment`.
#' @param label_col1 First label column name.
#' @param label_col2 Second label column name.
#' @param drop_na Logical; if TRUE, rows with NA in either label are dropped. Default TRUE.
#'
#' @return A data frame with columns `label`, `match_n`, `union_n`, and `agreement`.
#' @export
label_agreement_rates <- function(spe,
                                  label_col1,
                                  label_col2,
                                  drop_na = TRUE) {
  .assert_spe(spe)
  cd <- SummarizedExperiment::colData(spe)

  if (!label_col1 %in% colnames(cd)) {
    stop("label_col1 '", label_col1, "' not found in colData(spe).")
  }
  if (!label_col2 %in% colnames(cd)) {
    stop("label_col2 '", label_col2, "' not found in colData(spe).")
  }

  v1 <- cd[[label_col1]]
  v2 <- cd[[label_col2]]

  if (drop_na) {
    keep <- !(is.na(v1) | is.na(v2))
    v1 <- v1[keep]
    v2 <- v2[keep]
  }

  labels_all <- sort(unique(c(as.character(v1), as.character(v2))))
  v1 <- factor(v1, levels = labels_all)
  v2 <- factor(v2, levels = labels_all)

  agree_df <- data.frame(label = labels_all, stringsAsFactors = FALSE)
  agree_df$match_n <- vapply(labels_all, function(lbl) {
    sum(v1 == lbl & v2 == lbl, na.rm = TRUE)
  }, integer(1))
  agree_df$union_n <- vapply(labels_all, function(lbl) {
    sum(v1 == lbl | v2 == lbl, na.rm = TRUE)
  }, integer(1))
  agree_df$agreement <- ifelse(agree_df$union_n > 0,
                               agree_df$match_n / agree_df$union_n,
                               NA_real_)

  agree_df
}




#' Predict and fill unknown labels with a trained model
#'
#' Uses a trained `ranger` model to predict labels for all cells and fills
#' missing or `unknown_label` entries in `colData`.
#'
#' @param spe A `SpatialExperiment` or `SingleCellExperiment`.
#' @param model A fitted `ranger` model returned by [train_custom_rf()].
#' @param assay_name Assay to use as features (default "norm").
#' @param label_col Name of the existing label column (default "custom_label").
#' @param out_col Name of the output column with filled labels (default "soft_tree_label_filled").
#' @param pred_col Optional column to store raw model predictions (default "rf_pred").
#' @param unknown_label Label value to treat as missing (default "Unknown").
#'
#' @return The input `spe` with updated `colData` columns.
#' @export
predict_unknown_with_rf <- function(spe,
                                    model,
                                    assay_name = "norm",
                                    label_col = "custom_label",
                                    out_col = "soft_tree_label_filled",
                                    pred_col = "rf_pred",
                                    unknown_label = "Unknown") {
  .assert_spe(spe)
  if (!inherits(model, "ranger")) stop("model must be a ranger object.")
  if (!label_col %in% colnames(SummarizedExperiment::colData(spe))) {
    stop("label_col not found in colData(spe).")
  }
  if (!assay_name %in% SummarizedExperiment::assayNames(spe)) {
    stop(paste0("Assay '", assay_name, "' not found in spe."))
  }

  feat_mat <- SummarizedExperiment::assay(spe, assay_name)
  feature_df <- as.data.frame(t(feat_mat))

  # Align feature columns with the model's expected variables
  expected <- model$forest$independent.variable.names
  missing_cols <- setdiff(expected, names(feature_df))
  if (length(missing_cols) > 0) {
    stop("Missing feature columns required by model: ", paste(missing_cols, collapse = ", "))
  }
  feature_df <- feature_df[, expected, drop = FALSE]

  preds <- predict(model, feature_df)$predictions

  labels <- SummarizedExperiment::colData(spe)[[label_col]]
  filled <- labels
  replace_idx <- is.na(labels) | labels == unknown_label
  filled[replace_idx] <- as.character(preds)[replace_idx]

  SummarizedExperiment::colData(spe)[[pred_col]] <- preds
  SummarizedExperiment::colData(spe)[[out_col]] <- filled
  spe
}


#' Tabular and text summaries of RF performance
#'
#' Convenience helpers to expose overall metrics, cross-validation summaries,
#' and per-class statistics in both tabular and text-friendly formats.
#'
#' @param fit List returned by [train_custom_rf()].
#' @param digits Number of digits when formatting text output (default 3).
#'
#' @return `rf_metric_table()` returns a list of tables: `overall`,
#'   `cv_overall`, `cv_overall_summary`, `class`, and `class_raw`.
#'   `rf_metric_text()` returns a character vector of summary lines.
#' @export
rf_metric_table <- function(fit) {
  if (is.null(fit$metrics)) stop("fit must contain a metrics element.")

  metrics <- fit$metrics

  overall <- data.frame(
    metric = c("accuracy", "f1_macro"),
    value = c(metrics$accuracy, metrics$f1_macro),
    stringsAsFactors = FALSE
  )

  cv_overall_summary <- NULL
  if (!is.null(metrics$cv_overall) && nrow(metrics$cv_overall) > 0) {
    cv_overall_summary <- metrics$cv_overall |>
      dplyr::summarise(
        accuracy_mean = mean(accuracy, na.rm = TRUE),
        accuracy_sd = stats::sd(accuracy, na.rm = TRUE),
        f1_macro_mean = mean(f1_macro, na.rm = TRUE),
        f1_macro_sd = stats::sd(f1_macro, na.rm = TRUE)
      )
  }

  class_metrics <- class_metrics_from_fit(fit)

  class_summary <- class_metrics
  if ("fold" %in% colnames(class_metrics)) {
    class_summary <- class_metrics |>
      dplyr::group_by(class) |>
      dplyr::summarise(
        precision_mean = mean(precision, na.rm = TRUE),
        recall_mean = mean(recall, na.rm = TRUE),
        f1_mean = mean(f1, na.rm = TRUE),
        support_mean = mean(support, na.rm = TRUE),
        .groups = "drop"
      )
  }

  list(
    overall = overall,
    cv_overall = metrics$cv_overall,
    cv_overall_summary = cv_overall_summary,
    class = class_summary,
    class_raw = class_metrics
  )
}


#' @export
rf_metric_text <- function(fit, digits = 3) {
  tbl <- rf_metric_table(fit)
  fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))

  acc <- tbl$overall$value[tbl$overall$metric == "accuracy"]
  f1 <- tbl$overall$value[tbl$overall$metric == "f1_macro"]

  lines <- c(
    paste0("Overall accuracy: ", fmt(acc)),
    paste0("Overall F1 (macro): ", fmt(f1))
  )

  if (!is.null(tbl$cv_overall_summary)) {
    s <- tbl$cv_overall_summary
    lines <- c(
      lines,
      paste0("CV accuracy (mean ± sd): ", fmt(s$accuracy_mean), " ± ", fmt(s$accuracy_sd)),
      paste0("CV F1 macro (mean ± sd): ", fmt(s$f1_macro_mean), " ± ", fmt(s$f1_macro_sd))
    )
  }

  lines
}


#' Default quantile-based cutoff
#'
#' Computes a per-vector cutoff using the specified quantile.
#'
#' @param x Numeric vector of probabilities.
#' @param prob Quantile to use (default 0.98).
#'
#' @return Numeric scalar cutoff.
#' @export
prob_quantile_cutoff <- function(x, prob = 0.98) {
  x <- as.numeric(x)
  stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE)
}


#' Median + MAD cutoff helper
#'
#' Computes a cutoff as median(x) + `mad_mult` * MAD(x).
#'
#' @param x Numeric vector of probabilities.
#' @param mad_mult Multiplier on MAD (default 3).
#'
#' @return Numeric scalar cutoff.
#' @export
prob_mad_cutoff <- function(x, mad_mult = 3) {
  x <- as.numeric(x)
  med <- stats::median(x, na.rm = TRUE)
  mad <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)
  med + mad_mult * mad
}


#' Build a binary label matrix from probabilities
#'
#' Applies a per-cell-type cutoff function to each probability column to build
#' a logical matrix indicating which cell types pass their cutoff for each
#' cell.
#'
#' @param prob_mat Numeric probability matrix (cells x cell types).
#' @param cutoff_fn Function taking a numeric vector and returning a single
#'   numeric cutoff (default [prob_quantile_cutoff()]).
#'
#' @return Logical matrix with the same dimensions and dimnames as `prob_mat`.
#' @export
probability_label_matrix <- function(prob_mat, cutoff_fn = prob_quantile_cutoff) {
  if (!is.matrix(prob_mat)) stop("prob_mat must be a matrix.")
  if (nrow(prob_mat) == 0 || ncol(prob_mat) == 0) {
    stop("prob_mat must have at least one row and one column.")
  }
  if (!is.function(cutoff_fn)) stop("cutoff_fn must be a function.")

  cutoffs <- vapply(seq_len(ncol(prob_mat)), function(j) {
    val <- cutoff_fn(prob_mat[, j])
    if (length(val) != 1 || !is.numeric(val)) {
      stop("cutoff_fn must return a single numeric value per column.")
    }
    as.numeric(val)
  }, numeric(1))

  names(cutoffs) <- colnames(prob_mat)
  cutoffs[!is.finite(cutoffs)] <- NA_real_

  label_mat <- sweep(prob_mat, 2, cutoffs, FUN = ">=")
  mode(label_mat) <- "logical"
  label_mat
}


#' Apply cutoff-based labels to a result list
#'
#' Takes a `res` object (as returned by `run_soft_gating()` or
#' `run_tree_gating()`), builds a per-cell-type label matrix using a cutoff
#' function, stores it as `res$label_mat`, and writes a single-label assignment
#' into `colData(res$spe)`.
#'
#' A cell receives a label only if exactly one cell type passes its cutoff; if
#' multiple (or zero) cell types pass, the cell is labeled `unknown_label`.
#'
#' @param res List containing at least `spe` and `prob_mat`.
#' @param cutoff_fn Function passed to [probability_label_matrix()] to compute
#'   per-column cutoffs (default [prob_quantile_cutoff()]).
#' @param label_col Name of the column to add to `colData(res$spe)` (default
#'   "cutoff_label").
#' @param unknown_label Label used when zero or multiple cell types pass their
#'   cutoff (default "Unknown").
#'
#' @return The input `res` list with updated `spe` and a new `label_mat`
#'   element.
#' @export
apply_cutoff_labels <- function(res,
                                cutoff_fn = prob_quantile_cutoff,
                                label_col = "cutoff_label",
                                unknown_label = "Unknown") {
  if (is.null(res) || is.null(res$prob_mat)) stop("res$prob_mat is required.")
  if (is.null(res$spe)) stop("res$spe is required.")

  prob_mat <- res$prob_mat
  spe <- res$spe

  if (!is.matrix(prob_mat)) stop("res$prob_mat must be a matrix.")
  if (nrow(prob_mat) != ncol(spe)) {
    stop("nrow(res$prob_mat) must match ncol(res$spe).")
  }
  if (!is.null(rownames(prob_mat))) {
    if (!identical(rownames(prob_mat), colnames(spe))) {
      stop("rownames(res$prob_mat) must match colnames(res$spe) when provided.")
    }
  }

  label_mat <- probability_label_matrix(prob_mat, cutoff_fn = cutoff_fn)

  positive_counts <- rowSums(label_mat, na.rm = TRUE)
  labels <- rep(unknown_label, nrow(label_mat))

  single_idx <- which(positive_counts == 1)
  if (length(single_idx)) {
    single_labels <- apply(label_mat[single_idx, , drop = FALSE], 1, function(x) {
      colnames(label_mat)[which(x)[1]]
    })
    labels[single_idx] <- single_labels
  }

  SummarizedExperiment::colData(spe)[[label_col]] <- labels

  res$label_mat <- label_mat
  res$spe <- spe
  res
}



