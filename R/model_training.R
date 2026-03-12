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
                                      parallel = TRUE) {

  if (!is.null(seed)) set.seed(seed)
  .assert_spe(spe)

  num_threads <- if (parallel) max(1, parallel::detectCores(logical = FALSE)) else 1

  # 1. Feature Prep
  feat_mat_all <- SummarizedExperiment::assay(spe, assay_name)
  features_use <- if (length(features) == 1 && features == "all") rownames(feat_mat_all) else intersect(features, rownames(feat_mat_all))

  lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
  feature_df <- as.data.frame(t(feat_mat_all[features_use, , drop = FALSE]))

  # Filter only labelled cells for training/cleaning
  core_idx <- which(lab_vec != unknown_label & !is.na(lab_vec))
  core_df <- feature_df[core_idx, ]
  core_df$original_label <- factor(lab_vec[core_idx])

  # --- HELPER: F1 MACRO ---
  f1_macro <- function(truth, pred) {
    classes <- levels(truth)
    f1s <- sapply(classes, function(cls) {
      tp <- sum(truth == cls & pred == cls)
      fp <- sum(truth != cls & pred == cls)
      fn <- sum(truth == cls & pred != cls)
      prec <- if(tp + fp == 0) 0 else tp / (tp + fp)
      rec  <- if(tp + fn == 0) 0 else tp / (tp + fn)
      if(prec + rec == 0) 0 else 2 * (prec * rec) / (prec + rec)
    })
    return(f1s) # Returns vector for all classes
  }

  # --- STAGE 1: CONSENSUS CLEANING LOOP ---
  message(sprintf("Stage 1: Cleaning labels via %d repeats...", repeats))
  match_counts <- setNames(numeric(nrow(core_df)), rownames(core_df))

  for (r in seq_len(repeats)) {
    fold_assign <- integer(nrow(core_df))
    for (cls in levels(core_df$original_label)) {
      cls_idx <- which(core_df$original_label == cls)
      fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
    }
    for (k in seq_len(cv_folds)) {
      train_idx <- which(fold_assign != k); test_idx <- which(fold_assign == k)
      tmp <- ranger::ranger(original_label ~ ., data = core_df[train_idx, ], num.trees = 100, num.threads = num_threads)
      preds <- stats::predict(tmp, core_df[test_idx, ])$predictions
      match_counts[test_idx] <- match_counts[test_idx] + as.numeric(preds == core_df$original_label[test_idx])
    }
  }

  agreement_rate <- match_counts / repeats
  valid_names <- names(agreement_rate)[agreement_rate >= agreement_thresh]
  cleaned_df <- core_df[valid_names, ]

  # --- STAGE 2: EVALUATION LOOP (On Cleaned Data Only) ---
  message("Stage 2: Evaluating final model performance on cleaned cells...")
  eval_preds <- factor(rep(NA_character_, nrow(cleaned_df)), levels = levels(cleaned_df$original_label))
  eval_fold_assign <- integer(nrow(cleaned_df))
  for (cls in levels(cleaned_df$original_label)) {
    cls_idx <- which(cleaned_df$original_label == cls)
    eval_fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
  }

  overall_list <- list()
  class_metric_list <- list()

  for (k in seq_len(cv_folds)) {
    train_idx <- which(eval_fold_assign != k); test_idx <- which(eval_fold_assign == k)
    eval_mod <- ranger::ranger(original_label ~ ., data = cleaned_df[train_idx, ], num.trees = num.trees, num.threads = num_threads)
    pk <- stats::predict(eval_mod, cleaned_df[test_idx, ])$predictions
    pk <- factor(pk, levels = levels(cleaned_df$original_label))
    tk <- cleaned_df$original_label[test_idx]

    eval_preds[test_idx] <- pk

    # Calculate Overall Fold Metrics
    f1_vec <- f1_macro(tk, pk)
    overall_list[[k]] <- data.frame(fold = k, accuracy = mean(pk == tk), f1_macro = mean(f1_vec))

    # Calculate Per-Class Metrics for this fold
    classes <- levels(cleaned_df$original_label)
    class_metric_list[[k]] <- data.frame(
      class = classes,
      fold = k,
      precision = sapply(classes, function(c) {
        tp <- sum(tk==c & pk==c); fp <- sum(tk!=c & pk==c)
        if(tp+fp==0) 0 else tp/(tp+fp)
      }),
      recall = sapply(classes, function(c) {
        tp <- sum(tk==c & pk==c); fn <- sum(tk==c & pk!=c)
        if(tp+fn==0) 0 else tp/(tp+fn)
      }),
      f1 = f1_vec
    )
  }

  # --- FINAL MODEL TRAINING ---
  mtry_val <- if (is.null(mtry)) floor(sqrt(length(features_use))) else mtry
  final_model <- ranger::ranger(original_label ~ ., data = cleaned_df,
                                num.trees = num.trees, mtry = mtry_val,
                                probability = TRUE, num.threads = num_threads)

  # Update SPE colData
  new_labels <- lab_vec
  new_labels[colnames(spe) %in% names(agreement_rate[agreement_rate < agreement_thresh])] <- unknown_label
  SummarizedExperiment::colData(spe)$cleaned_core_label <- new_labels

  # --- RETURN LIST (Matches your old output structure) ---
  return(list(
    spe = spe,
    model = final_model,
    agreement_rates = agreement_rate,
    features_used = features_use,
    metrics = list(
      confusion_matrix = table(truth = cleaned_df$original_label, pred = eval_preds),
      accuracy = mean(cleaned_df$original_label == eval_preds),
      f1_macro = mean(f1_macro(cleaned_df$original_label, eval_preds)),
      cv_overall = dplyr::bind_rows(overall_list),
      cv_class_metrics = dplyr::bind_rows(class_metric_list)
    ),
    test_pred = eval_preds,
    test_truth = cleaned_df$original_label
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

  return(spe)
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



