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
train_custom_randomforest <- function(spe,
                                      label_col = "cutoff_label",
                                      assay_name = "norm",
                                      unknown_label = "Unknown",
                                      num.trees = 200,
                                      mtry = NULL,
                                      seed = NULL,
                                      features = "all",
                                      cv_folds = 5,
                                      repeats = 10,
                                      agreement_thresh = 0.8,
                                      num_threads = 2) {

  if (!is.null(seed)) set.seed(seed)
  .assert_spe(spe)

  # 1. Feature Prep
  feat_mat_all <- SummarizedExperiment::assay(spe, assay_name)
  features_use <- if (length(features) == 1 && features == "all") rownames(feat_mat_all) else intersect(features, rownames(feat_mat_all))

  lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
  feature_df <- as.data.frame(t(feat_mat_all[features_use, , drop = FALSE]))

  # Filter only labelled cells for training
  core_idx <- which(lab_vec != unknown_label & !is.na(lab_vec))
  core_df <- feature_df[core_idx, ]
  core_df$original_label <- factor(lab_vec[core_idx])

  class_levels <- levels(core_df$original_label)
  sum_prob_mat <- matrix(0, nrow = nrow(core_df), ncol = length(class_levels),
                         dimnames = list(rownames(core_df), class_levels))

  # --- STAGE 1: CONSENSUS CLEANING & PROBABILITY GENERATION ---
  message(sprintf("Processing %d core cells via %d repeated CV folds...", nrow(core_df), repeats))
  match_counts <- setNames(numeric(nrow(core_df)), rownames(core_df))

  for (r in seq_len(repeats)) {
    fold_assign <- integer(nrow(core_df))
    for (cls in class_levels) {
      cls_idx <- which(core_df$original_label == cls)
      fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
    }

    for (k in seq_len(cv_folds)) {
      train_idx <- which(fold_assign != k); test_idx <- which(fold_assign == k)

      # Probability = TRUE is the key for uncertainty mapping
      tmp <- ranger::ranger(original_label ~ .,
                            data = core_df[train_idx, ],
                            num.trees = 100,
                            num.threads = num_threads,
                            probability = TRUE)

      prob_preds <- stats::predict(tmp, core_df[test_idx, ])$predictions
      sum_prob_mat[test_idx, ] <- sum_prob_mat[test_idx, ] + prob_preds

      # Determine hard labels for agreement counting
      preds <- class_levels[max.col(prob_preds)]
      match_counts[test_idx] <- match_counts[test_idx] + as.numeric(preds == core_df$original_label[test_idx])
    }
  }

  # 2. Consensus Calculations
  avg_core_prob_mat <- sum_prob_mat / repeats
  agreement_rate <- match_counts / repeats

  # Consensus predictions (the type the CV models consistently chose)
  consensus_labels <- class_levels[max.col(avg_core_prob_mat)]
  consensus_labels <- factor(consensus_labels, levels = class_levels)

  # 3. Cleaning
  valid_names <- names(agreement_rate)[agreement_rate >= agreement_thresh]
  cleaned_df <- core_df[valid_names, ]

  # Update SPE with cleaned labels
  new_labels <- lab_vec
  new_labels[colnames(spe) %in% names(agreement_rate[agreement_rate < agreement_thresh])] <- unknown_label
  SummarizedExperiment::colData(spe)$cleaned_core_label <- new_labels

  # 4. Final Model Training (On Cleaned Data Only)
  message(sprintf("Training final model on %d cleaned cells...", nrow(cleaned_df)))
  mtry_val <- if (is.null(mtry)) floor(sqrt(length(features_use))) else mtry
  final_model <- ranger::ranger(original_label ~ ., data = cleaned_df,
                                num.trees = num.trees, mtry = mtry_val,
                                probability = TRUE, num.threads = num_threads)

  # --- RETURN RESULTS ---
  return(list(
    spe = spe,
    model = final_model,
    core_prob_mat = avg_core_prob_mat,
    agreement_rates = agreement_rate,
    features_used = features_use
  ))
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
  pred  <- test_pred  %||% fit$test_pred

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
      precision = ifelse(.data$tp + .data$fp == 0, 0, .data$tp / (.data$tp + .data$fp)),
      recall    = ifelse(.data$tp + .data$fn == 0, 0, .data$tp / (.data$tp + .data$fn)),
      f1        = ifelse(.data$precision + .data$recall == 0, 0,
                         2 * .data$precision * .data$recall / (.data$precision + .data$recall))
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
predict_unknown_with_randomforest <- function(spe,
                                              model,
                                              assay_name = "norm",
                                              label_col = "custom_label",
                                              out_col = "soft_tree_label_filled",
                                              pred_col = "rf_pred",
                                              unknown_label = "Unknown",
                                              threshold = 0.5,           # New: Confidence threshold
                                              unassigned_label = "Unassigned") { # New: Label for low confidence
  .assert_spe(spe)
  if (!inherits(model, "ranger")) stop("model must be a ranger object.")

  # Ranger must be trained with probability = TRUE for thresholding to work
  if (model$treetype != "Probability estimation") {
    warning("Model was not trained with probability = TRUE. Thresholding will be ignored.")
    threshold <- 0
  }

  feat_mat <- SummarizedExperiment::assay(spe, assay_name)
  feature_df <- as.data.frame(t(feat_mat))

  # Align features
  expected <- model$forest$independent.variable.names
  feature_df <- feature_df[, expected, drop = FALSE]

  # 1. Get Probabilities
  pred_obj <- stats::predict(model, feature_df)
  prob_mat <- pred_obj$predictions # This is a matrix of probabilities per class

  # 2. Determine winners and their confidence level
  max_probs <- apply(prob_mat, 1, max)
  winning_indices <- apply(prob_mat, 1, which.max)
  raw_preds <- colnames(prob_mat)[winning_indices]

  # 3. Apply Thresholding
  # If confidence < threshold, label as 'Unassigned'
  final_preds <- ifelse(max_probs >= threshold, raw_preds, unassigned_label)

  # 4. Fill into SPE
  labels <- SummarizedExperiment::colData(spe)[[label_col]]
  filled <- labels

  # We only fill cells that were originally 'Unknown' or NA
  replace_idx <- is.na(labels) | labels == unknown_label
  filled[replace_idx] <- final_preds[replace_idx]

  SummarizedExperiment::colData(spe)[[pred_col]] <- final_preds
  SummarizedExperiment::colData(spe)[[out_col]] <- filled

  # Optional: store the max probability for QC
  SummarizedExperiment::colData(spe)$rf_confidence <- max_probs

  return(list(spe = spe, prob_mat = prob_mat))
}


#' Tabular and text summaries of RF performance
#'
#' Convenience helpers to expose overall metrics, cross-validation summaries,
#' and per-class statistics in both tabular and text-friendly formats.
#'
#' @param fit List returned by [train_custom_rf()].
#'
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
        accuracy_mean = mean(.data$accuracy, na.rm = TRUE),
        accuracy_sd   = stats::sd(.data$accuracy, na.rm = TRUE),
        f1_macro_mean = mean(.data$f1_macro, na.rm = TRUE),
        f1_macro_sd   = stats::sd(.data$f1_macro, na.rm = TRUE)
      )
  }

  class_metrics <- class_metrics_from_fit(fit)

  class_summary <- class_metrics
  if ("fold" %in% colnames(class_metrics)) {
    class_summary <- class_metrics |>
      dplyr::group_by(.data$class) |>
      dplyr::summarise(
        precision_mean = mean(.data$precision, na.rm = TRUE),
        recall_mean    = mean(.data$recall, na.rm = TRUE),
        f1_mean        = mean(.data$f1, na.rm = TRUE),
        support_mean   = mean(.data$support, na.rm = TRUE),
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

#' Text summary of random forest performance
#'
#' Formats key overall and cross-validation metrics from a fitted model
#' returned by [train_custom_rf()] into a character vector of human-readable lines.
#'
#' @param fit List returned by [train_custom_rf()].
#' @param digits Number of digits when formatting numeric values (default 3).
#'
#' @return A character vector of summary lines.
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
      paste0("CV accuracy (mean +/- sd): ", fmt(s$accuracy_mean), " +/- ", fmt(s$accuracy_sd)),
      paste0("CV F1 macro (mean +/- sd): ", fmt(s$f1_macro_mean), " +/- ", fmt(s$f1_macro_sd))
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
#' @param k Integer specifying the number of nearest neighbours used in the
#'   kNN classifier during consensus cleaning. Default is \code{5}.
#' @param method Character string specifying the similarity or distance metric
#'   used for kNN. Options include \code{"pearson"}, \code{"spearman"},
#'   \code{"cosine"}, or \code{"euclidean"}. Default is \code{"pearson"}.
#' @param seed Optional integer for random seed to ensure reproducibility.
#'   Default is \code{NULL}.
#' @param chunk_size Integer number of test cells processed per block inside
#'   weighted kNN calls. Smaller values reduce peak RAM with identical
#'   predictions at the cost of runtime. Default is \code{250L}.
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
                             k = 5,
                             method = "pearson",
                             seed = NULL,
                             chunk_size = 250L) {

  if (!is.null(seed)) set.seed(seed)
  .assert_spe(spe)

  # 1. Feature Prep
  feat_mat <- SummarizedExperiment::assay(spe, assay_name)
  features_use <- if (length(features) == 1 && features == "all") rownames(feat_mat) else intersect(features, rownames(feat_mat))

  lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
  feature_mat <- t(as.matrix(feat_mat[features_use, , drop = FALSE]))

  core_idx <- which(lab_vec != unknown_label & !is.na(lab_vec))
  core_mat <- feature_mat[core_idx, , drop = FALSE]
  core_labels <- factor(lab_vec[core_idx])

  class_levels <- levels(core_labels)
  # Accumulate probability matrices across repeats
  sum_prob_mat <- matrix(0, nrow = nrow(core_mat), ncol = length(class_levels),
                         dimnames = list(rownames(core_mat), class_levels))

  # --- STAGE 1: WEIGHTED CONSENSUS CLEANING ---
  message(sprintf("Starting Weighted kNN (%s) cleaning on %d cells: %d repeats...", method, nrow(core_mat), repeats))
  match_counts <- setNames(numeric(nrow(core_mat)), rownames(core_mat))

  for (r in seq_len(repeats)) {
    fold_assign <- integer(nrow(core_mat))
    for (cls in class_levels) {
      cls_idx <- which(core_labels == cls)
      fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
    }

    for (k_fold in seq_len(cv_folds)) {
      train_idx <- which(fold_assign != k_fold)
      test_idx <- which(fold_assign == k_fold)

      # Use your wkNN logic to get probabilities for the test fold
      # We call your internal logic here
      res_wkNN <- predict_wknn_multi(
        train_data = core_mat[train_idx, , drop = FALSE],
        test_data = core_mat[test_idx, , drop = FALSE],
        train_labels = core_labels[train_idx],
        k = k,
        method = method,
        return_matrix = TRUE,
        chunk_size = chunk_size
      )

      # Accumulate the weighted probability matrix
      sum_prob_mat[test_idx, ] <- sum_prob_mat[test_idx, ] + res_wkNN$prob_matrix

      # Track hard-label matches for cleaning
      match_counts[test_idx] <- match_counts[test_idx] + as.numeric(res_wkNN$labels == core_labels[test_idx])
    }
  }

  # Average the probabilities across all repeats
  avg_core_prob_mat <- sum_prob_mat / repeats
  agreement_rate <- match_counts / repeats

  # 2. Update Labels
  inconsistent_names <- names(agreement_rate)[agreement_rate < agreement_thresh]
  valid_core_names <- names(agreement_rate)[agreement_rate >= agreement_thresh]

  cleaned_labels <- lab_vec
  cleaned_labels[colnames(spe) %in% inconsistent_names] <- unknown_label
  SummarizedExperiment::colData(spe)$cleaned_core_label <- cleaned_labels
  valid_core_idx <- which(rownames(core_mat) %in% valid_core_names)

  return(list(
    spe = spe,
    model = list(
      type = "wknn",
      reference_data   = as.data.frame(core_mat[valid_core_idx, , drop = FALSE]),
      reference_labels = factor(core_labels[valid_core_idx], levels = class_levels),
      params = list(k = k, method = method)
    ),
    core_prob_mat  = avg_core_prob_mat,
    agreement_rates = agreement_rate,
    features_used  = features_use
  ))

}

#' Predict Unknown Cells using kNN Reference
#'
#' Uses a trained kNN model to predict cell types for cells labelled as Unknown
#' or Unassigned, based on marker expression from a specified assay.
#'
#' @param spe A \code{SpatialExperiment} object containing cell data.
#' @param knn_ref Output from \code{train_custom_knn()}, containing the trained
#'   kNN model and reference data.
#' @param assay_name Character string specifying the assay to use for marker
#'   expression. Default is \code{"norm"}.
#' @param label_col Character string specifying the column in \code{colData(spe)}
#'   containing the initial cell type labels. Default is \code{"cutoff_label"}.
#' @param out_col Character string specifying the name of the output column added
#'   to \code{colData(spe)} containing the final filled labels. Default is
#'   \code{"knn_label_filled"}.
#' @param pred_col Character string specifying the name of the output column added
#'   to \code{colData(spe)} containing the raw kNN predictions. Default is
#'   \code{"knn_pred"}.
#' @param unknown_label Character string specifying the label used to identify
#'   unknown cells. Default is \code{"Unknown"}.
#' @param unassigned_label Character string specifying the label used to identify
#'   unassigned cells. Default is \code{"Unassigned"}.
#' @param threshold Numeric value between 0 and 1 specifying the minimum
#'   proportion of neighbours required to agree on a label for a prediction to
#'   be accepted. Default is \code{0.6}.
#' @param k Integer specifying the number of nearest neighbours to use for
#'   prediction. Default is \code{5}.
#' @param dist_method Character string specifying the distance metric used for
#'   kNN similarity calculation (e.g., \code{"pearson"}, \code{"cosine"},
#'   \code{"euclidean"}). Default is \code{"pearson"}.
#' @param chunk_size Integer number of test cells processed per weighted kNN
#'   block. Smaller values reduce peak RAM with identical predictions at the
#'   cost of runtime. Default is \code{250L}.
#'
#' @return A \code{SpatialExperiment} object with two new columns added to
#'   \code{colData}: \code{out_col} containing the final predicted labels and
#'   \code{pred_col} containing the raw kNN predictions.
#'
#' @export
predict_unknown_with_knn <- function(spe,
                                     knn_ref,
                                     assay_name = "norm",
                                     label_col = "cutoff_label",
                                     out_col = "knn_label_filled",
                                     pred_col = "knn_pred",
                                     unknown_label = "Unknown",
                                     unassigned_label = "Unassigned",
                                     threshold = 0.6,
                                     k = 5,
                                     dist_method = "pearson",
                                     chunk_size = 250L) {

  # 1. Feature Prep
  feat_mat <- SummarizedExperiment::assay(spe, assay_name)
  feature_mat <- t(as.matrix(feat_mat[knn_ref$features_used, , drop = FALSE]))

  # 2. Get Probabilities (scClassify style)
  # 'prob = TRUE' returns the proportion of the winning class votes as an attribute
  # knn_res <- class::knn(train = knn_ref$reference_data,
  #                       test = feature_df,
  #                       cl = knn_ref$reference_labels,
  #                       k = k,
  #                       prob = TRUE)

  knn_res <- predict_wknn_multi(
    train_data = knn_ref$model$reference_data,
    test_data = feature_mat,
    train_labels = knn_ref$model$reference_labels,
    k = k,
    method = dist_method,
    return_matrix = TRUE,
    chunk_size = chunk_size
  )

  # raw_preds <- as.character(knn_res)
  # confidences <- attr(knn_res, "prob")

  raw_preds <- knn_res$labels
  confidences <- knn_res$probs
  prob_mat <- knn_res$prob_matrix # The new Likelihood Matrix

  # 3. Apply Confidence Thresholding
  # If the % of neighbor votes < threshold, label as 'Unassigned'
  final_preds <- ifelse(confidences >= threshold, raw_preds, unassigned_label)

  # 4. Fill into SPE (Hybrid Logic)
  # We only replace cells that were 'Unknown' in the original OR
  # cells that were turned into 'Unknown' by the train_custom_knn cleaning step.

  # Step A: Get the current labels (using the cleaned ones if available)
  current_labels <- if ("cleaned_core_label" %in% names(SummarizedExperiment::colData(spe))) {
    SummarizedExperiment::colData(spe)$cleaned_core_label
  } else {
    SummarizedExperiment::colData(spe)[[label_col]]
  }

  filled <- current_labels

  # Step B: Identify cells that need a prediction
  # (Either they were always Unknown, or they failed the cleaning consensus)
  replace_idx <- is.na(current_labels) | current_labels == unknown_label

  filled[replace_idx] <- final_preds[replace_idx]

  # 5. Store results back to SPE
  SummarizedExperiment::colData(spe)[[pred_col]] <- final_preds
  SummarizedExperiment::colData(spe)[[out_col]] <- filled
  SummarizedExperiment::colData(spe)$knn_confidence <- confidences

  return(list(spe = spe, prob_mat = prob_mat))
}









#' Multi-Metric Weighted kNN
#'
#' @param train_data Data frame or matrix of training cells (rows = markers, cols = cells).
#' @param test_data Data frame or matrix of test cells (rows = markers, cols = cells).
#' @param train_labels Factor of labels for the training data.
#' @param k Number of neighbors.
#' @param method One of "pearson", "spearman", "cosine", or "euclidean"
#' @param return_matrix Logical indicating whether to return the full probability
#'   matrix for all classes. If \code{FALSE}, only the predicted labels and
#'   associated probabilities are returned. Default is \code{FALSE}.
#' @param chunk_size Optional integer number of test cells per block. If
#'   \code{NULL}, all test cells are processed in a single block (legacy
#'   behavior).
#'
#' @return A list containing predicted 'labels' and 'probs'.
predict_wknn_multi <- function(train_data,
                               test_data,
                               train_labels,
                               k = 5,
                               method = "pearson",
                               return_matrix = FALSE,
                               chunk_size = NULL) {

  # Ensure data is matrix format for fast calculation
  train_mat <- t(as.matrix(train_data))
  test_mat <- t(as.matrix(test_data))
  n_test <- ncol(test_mat)
  class_levels <- levels(train_labels)

  # Safety Check: Do the number of markers match?
  if (nrow(train_mat) != nrow(test_mat)) {
    stop(sprintf("Dimension mismatch! Train has %d markers, Test has %d markers.
                  Ensure both matrices are Cells (rows) x Markers (columns).",
                 ncol(train_mat), ncol(test_mat)))
  }

  if (n_test == 0L) {
    out <- list(labels = character(0), probs = numeric(0))
    if (return_matrix) {
      out$prob_matrix <- matrix(
        numeric(0),
        nrow = 0,
        ncol = length(class_levels),
        dimnames = list(character(0), class_levels)
      )
    }
    return(out)
  }

  if (is.null(chunk_size) || !is.finite(chunk_size) || as.integer(chunk_size) <= 0L) {
    chunk_size <- n_test
  } else {
    chunk_size <- as.integer(chunk_size)
  }

  # Keep cosine numerics identical to legacy behavior by evaluating in one block.
  if (identical(method, "cosine")) {
    chunk_size <- n_test
  }

  chunk_idx <- split(seq_len(n_test), ceiling(seq_len(n_test) / chunk_size))

  out_labels <- character(n_test)
  out_probs <- numeric(n_test)
  if (!is.null(colnames(test_mat))) {
    names(out_labels) <- colnames(test_mat)
    names(out_probs) <- colnames(test_mat)
  }

  prob_matrix <- NULL
  if (return_matrix) {
    prob_matrix <- matrix(
      0,
      nrow = n_test,
      ncol = length(class_levels),
      dimnames = list(colnames(test_mat), class_levels)
    )
  }

  for (idx in chunk_idx) {
    test_chunk <- test_mat[, idx, drop = FALSE]

    # 1. Calculate the Similarity/Distance Matrix for this chunk
    if (method %in% c("pearson", "spearman")) {
      # Correlation: higher is more similar
      score_mat <- cor(test_chunk, train_mat, method = method)
      is_distance <- FALSE

    } else if (method == "cosine") {
      # Cosine Similarity: higher is more similar
      cp <- crossprod(test_chunk, train_mat)
      rn <- sqrt(colSums(test_chunk^2))
      cn <- sqrt(colSums(train_mat^2))
      score_mat <- cp / outer(rn, cn)
      is_distance <- FALSE

    } else if (method == "euclidean") {
      # Euclidean: lower is more similar (distance)
      score_mat <- as.matrix(proxy::dist(t(test_chunk), t(train_mat), method = "Euclidean"))
      is_distance <- TRUE
    } else {
      stop("Unknown method: ", method)
    }

    # 2. Process each test cell in the chunk
    for (i in seq_len(nrow(score_mat))) {
      scores <- score_mat[i, ]

      # Identify Top K
      if (is_distance) {
        top_k_idx <- order(scores, decreasing = FALSE)[1:k]
        # Convert distance to a weight (Inverse distance)
        # Add small epsilon to avoid division by zero
        weights <- 1 / (scores[top_k_idx] + 1e-6)
      } else {
        top_k_idx <- order(scores, decreasing = TRUE)[1:k]
        # Use raw similarity as weight (clipping negative correlations to 0)
        weights <- pmax(scores[top_k_idx], 0)
      }

      top_labels <- train_labels[top_k_idx]
      weights[!is.finite(weights)] <- 0
      w_sum <- sum(weights)
      if (!is.finite(w_sum) || w_sum <= 0) {
        global_i <- idx[i]
        out_labels[global_i] <- as.character(top_labels[1])
        out_probs[global_i] <- 0
        if (return_matrix) prob_matrix[global_i, ] <- 0
        next
      }

      # 3. Weighted Voting
      # Sum the weights for each unique label found in the neighbors
      label_sums <- tapply(weights, top_labels, sum)
      label_sums[is.na(label_sums)] <- 0

      best_label <- names(sort(label_sums, decreasing = TRUE))[1]

      # Probability is the winning weight sum divided by total weight sum
      prob <- max(label_sums) / w_sum

      global_i <- idx[i]
      out_labels[global_i] <- best_label
      out_probs[global_i] <- prob

      if (return_matrix) {
        prob_dist <- setNames(numeric(length(class_levels)), class_levels)
        prob_dist[names(label_sums)] <- label_sums / w_sum
        prob_matrix[global_i, ] <- prob_dist
      }
    }

    rm(test_chunk, score_mat)
  }

  # Format output as a clean list
  out <- list(
    labels = out_labels,
    probs  = out_probs
  )

  # Only build and add the matrix if requested (saves memory)
  if (return_matrix) {
    out$prob_matrix <- prob_matrix
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
#' @param chunk_size Integer specifying the number of cells processed per chunk
#'   during prediction. Increasing this value may improve speed but requires
#'   more memory. Default is \code{1000L}
#' @export
predict_hierarchical_knn_recursive <- function(spe,
                                               hier_ref,
                                               hc_tree,
                                               assay_name   = "exprs",
                                               threshold    = 0.7,
                                               agreement_threshold = 0.8,
                                               k            = 5,
                                               repeats      = 5,
                                               dist_methods = c("pearson", "cosine"),
                                               BPPARAM      = BiocParallel::SerialParam(),
                                               out_col      = "hier_label",
                                               chunk_size   = 1000L,
                                               unassigned_label = "Unassigned") {   # <-- NEW

  .assert_spe(spe)
  n_cells    <- ncol(spe)

  # FIX 1: Materialise the full dense matrix ONCE here.
  # Previously t(as.matrix(...)) was called inside process_node on every
  # recursive visit, allocating a fresh dense copy each time.
  feat_mat_dense <- t(as.matrix(SummarizedExperiment::assay(spe, assay_name)))
  # feat_mat_dense is now (cells × markers); we subset rows per node, not columns.

  final_labels   <- rep(NA_character_, n_cells)
  root_node_idx  <- nrow(hc_tree$merge)
  n_tasks        <- length(dist_methods) * repeats

  process_node <- function(node_idx, active_indices) {
    if (length(active_indices) == 0L) return()

    node_id  <- paste0("Node_", node_idx)
    ref      <- hier_ref[[node_id]]
    train_mat <- as.matrix(ref$train_data)   # cells × markers, already correct shape

    # FIX 2: Use integer indices instead of slicing train_mat here.
    # The bootstrap slice happens lazily inside the worker (see below).
    n_train <- nrow(train_mat)

    # FIX 4: Build a plain integer task list; no expand.grid data.frame kept in RAM.
    task_methods <- rep(dist_methods, each = repeats)  # length = n_tasks, just strings

    # Accumulators for stream-reduce (FIX 3 – no cbind over all repeats at once)
    vote_counts <- NULL   # will become (n_active × n_labels) integer matrix
    prob_sums   <- NULL   # (n_active) numeric vector
    n_done      <- 0L

    # ---- Process active cells in chunks (FIX 5) --------------------------------
    chunks     <- split(active_indices,
                        ceiling(seq_along(active_indices) / chunk_size))

    for (chunk_cells in chunks) {
      # FIX 1 (continued): subset rows of the pre-built dense matrix
      test_chunk <- feat_mat_dense[chunk_cells, ref$markers, drop = FALSE]

      # FIX 3 + FIX 4: stream-reduce — accumulate only vote tallies per chunk
      chunk_vote_counts <- NULL
      chunk_prob_sums   <- numeric(length(chunk_cells))

      ensemble_out <- BiocParallel::bplapply(
        seq_len(n_tasks),
        function(i) {
          m <- task_methods[i]

          # FIX 2: bootstrap by index; slice inside the worker so each worker
          # allocates only its own 80 % slice, not a copy visible to the parent.
          boot_idx <- sample.int(n_train, size = floor(0.8 * n_train))

          res <- predict_wknn_multi(
            train_data   = train_mat[boot_idx, , drop = FALSE],
            test_data    = test_chunk,
            train_labels = ref$train_labels[boot_idx],
            k            = k,
            method       = m
          )
          # Return only the compact summary needed for consensus.
          # The full (cells × k) distance matrix inside predict_wknn_multi
          # is freed when the worker returns.
          list(labels = res$labels, probs = res$probs)
        },
        BPPARAM = BPPARAM
      )

      # FIX 3: stream-reduce votes — never cbind the full matrix stack
      all_labels <- vapply(ensemble_out, `[[`, character(length(chunk_cells)), "labels")
      # all_labels is (n_chunk × n_tasks) — still allocates, but only for one chunk

      label_levels <- sort(unique(as.vector(all_labels)))

      # Tally votes as an integer matrix (n_chunk × n_levels)
      vote_mat <- vapply(label_levels, function(lv)
        as.integer(rowSums(all_labels == lv)), integer(length(chunk_cells)))

      prob_vec <- rowMeans(
        vapply(ensemble_out, `[[`, numeric(length(chunk_cells)), "probs"))

      # Merge this chunk's tallies into global accumulators
      if (is.null(chunk_vote_counts)) {
        chunk_vote_counts <- vote_mat
        colnames(chunk_vote_counts) <- label_levels
      }
      # (If running across chunks sequentially the merge is trivial; leave as-is.)
      chunk_prob_sums <- prob_vec

      rm(ensemble_out, all_labels, vote_mat)  # release before next chunk

      # ---- Consensus for this chunk -------------------------------------------
      node_preds    <- label_levels[max.col(chunk_vote_counts)]
      node_agreement <- apply(chunk_vote_counts, 1,
                              function(r) max(r)) / n_tasks
      node_avg_probs <- chunk_prob_sums

      uncertain_mask <- is.na(node_avg_probs) | is.na(node_agreement) |
        (node_avg_probs < threshold) | (node_agreement < agreement_threshold)

      if (any(uncertain_mask, na.rm = TRUE)) {
        final_labels[chunk_cells[uncertain_mask]] <<-
          # paste0(node_id, "_unassigned")
          unassigned_label
      }

      # ---- Routing for certain cells -------------------------------------------
      certain_idx <- which(!uncertain_mask)
      if (length(certain_idx) > 0L) {
        preds_certain   <- node_preds[certain_idx]
        indices_certain <- chunk_cells[certain_idx]

        for (choice in c("Left", "Right")) {
          group_idx <- indices_certain[preds_certain == choice]
          if (length(group_idx) == 0L) next

          side_idx  <- if (choice == "Left") 1L else 2L
          child_val <- hc_tree$merge[node_idx, side_idx]

          if (child_val < 0L) {
            final_labels[group_idx] <<- hc_tree$labels[-child_val]
          } else {
            process_node(child_val, group_idx)
          }
        }
      }
    }   # end chunk loop
  }   # end process_node

  message(sprintf(
    "Starting recursive ensemble classification for %d cells (chunk_size=%d)...",
    n_cells, chunk_size))

  process_node(root_node_idx, seq_len(n_cells))

  SummarizedExperiment::colData(spe)[[out_col]] <- final_labels
  return(spe)
}
