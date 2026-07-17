## deep_learning.R
##
## Deep-learning counterparts to train_custom_randomforest() / train_custom_knn()
## and predict_unknown_with_randomforest() / predict_unknown_with_knn().
##
## Requires the 'torch' package (Suggests, not Imports) -- all functions here
## check for it explicitly so the rest of the package works without it.


#' @keywords internal
.require_torch <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(
      "The 'torch' package is required for deep-learning functions.\n",
      "Install it with: install.packages('torch'); torch::install_torch()",
      call. = FALSE
    )
  }
}

#' Fit centering/scaling stats on a numeric matrix (cells x markers)
#' @keywords internal
.dl_scale_fit <- function(mat) {
  ctr <- colMeans(mat, na.rm = TRUE)
  scl <- apply(mat, 2, stats::sd)
  scl[!is.finite(scl) | scl == 0] <- 1  # avoid divide-by-zero for constant markers
  list(center = ctr, scale = scl)
}

#' @keywords internal
.dl_scale_apply <- function(mat, center, scale) {
  sweep(sweep(mat, 2, center, "-"), 2, scale, "/")
}


#' Define a small feedforward classifier: trunk of hidden layers + linear head.
#' Returns raw logits (not softmax) -- torch's cross-entropy loss expects logits.
#' @keywords internal
.build_dl_net <- function(input_dim, hidden_dims, n_classes, dropout) {
  torch::nn_module(
    initialize = function() {
      dims <- c(input_dim, hidden_dims)
      layers <- list()
      for (i in seq_along(hidden_dims)) {
        layers[[length(layers) + 1]] <- torch::nn_linear(dims[i], dims[i + 1])
        layers[[length(layers) + 1]] <- torch::nn_relu()
        layers[[length(layers) + 1]] <- torch::nn_dropout(p = dropout)
      }
      self$trunk <- do.call(torch::nn_sequential, layers)
      self$head  <- torch::nn_linear(dims[length(dims)], n_classes)
    },
    forward = function(x) {
      self$head(self$trunk(x))
    }
  )()
}

#' Simple rank-based (Mann-Whitney) binary AUC, no external dependency.
#' @keywords internal
.binary_auc <- function(scores, labels_binary) {
  pos <- scores[labels_binary == 1]
  neg <- scores[labels_binary == 0]
  n_pos <- length(pos); n_neg <- length(neg)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  sum_rank_pos <- sum(r[seq_len(n_pos)])
  (sum_rank_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}


#' Macro-averaged one-vs-rest AUC across classes (probs: cells x classes, 1-indexed y_int).
#' @keywords internal
.multiclass_auc <- function(prob_mat, y_int, n_classes) {
  aucs <- vapply(seq_len(n_classes), function(c) {
    .binary_auc(prob_mat[, c], as.integer(y_int == c))
  }, numeric(1))
  mean(aucs, na.rm = TRUE)
}


#' One eval-mode forward pass: returns loss, accuracy, macro OVR AUC for (x, y).
#' @keywords internal
.dl_eval_pass <- function(net, loss_fn, x, y_tensor, y_int_r, n_classes) {
  out <- torch::with_no_grad({ net(x) })
  loss <- as.numeric(loss_fn(out, y_tensor)$item())
  probs <- as.matrix(torch::as_array(torch::nnf_softmax(out, dim = 2)))
  preds <- max.col(probs)
  acc <- mean(preds == y_int_r)
  auc <- .multiclass_auc(probs, y_int_r, n_classes)
  list(loss = loss, acc = acc, auc = auc)
}

#'
#' Train one network from scratch on (x_mat, y_int), holding out a stratified
#' validation split each call. Tracks validation loss per epoch, keeps a copy
#' of the weights from the best epoch, restores them before returning, and
#' optionally stops early if validation loss hasn't improved for `patience`
#' epochs. y_int must be 1-indexed (1..n_classes) -- torch-for-R's cross-
#' entropy/nll losses require class targets starting at 1, not 0, despite the
#' underlying Python-ported docs saying [0, C-1].
#'
#' Returns a list: net (best-epoch weights, in eval() mode), history
#' (per-epoch train/valid loss+acc+auc data.frame), best_epoch, best_val_loss.
#' @keywords internal
.dl_train_loop <- function(x_mat, y_int, n_classes, hidden_dims, dropout,
                           epochs, lr, weight_decay, batch_size, class_weights = NULL,
                           val_frac = 0.15, patience = 10, verbose = FALSE) {
  net <- .build_dl_net(ncol(x_mat), hidden_dims, n_classes, dropout)
  opt <- torch::optim_adam(net$parameters, lr = lr, weight_decay = weight_decay)

  loss_fn <- if (!is.null(class_weights)) {
    torch::nn_cross_entropy_loss(weight = torch::torch_tensor(class_weights, dtype = torch::torch_float()))
  } else {
    torch::nn_cross_entropy_loss()
  }

  n <- nrow(x_mat)

  # Stratified train/validation split so every class is represented in both,
  # even for the smallest classes (at least 1 cell held out per class).
  val_idx <- integer(0)
  for (cls in seq_len(n_classes)) {
    cls_idx <- which(y_int == cls)
    if (length(cls_idx) == 0) next
    n_val_cls <- max(1, round(length(cls_idx) * val_frac))
    n_val_cls <- min(n_val_cls, length(cls_idx) - 1)  # keep at least 1 for training
    if (n_val_cls < 1) next
    val_idx <- c(val_idx, sample(cls_idx, n_val_cls))
  }
  train_idx <- setdiff(seq_len(n), val_idx)

  x_t <- torch::torch_tensor(as.matrix(x_mat), dtype = torch::torch_float())
  y_t <- torch::torch_tensor(as.integer(y_int), dtype = torch::torch_long())

  x_train <- x_t[train_idx, , drop = FALSE]
  y_train <- y_t[train_idx]
  x_val   <- x_t[val_idx, , drop = FALSE]
  y_val   <- y_t[val_idx]

  y_train_int_r <- y_int[train_idx]  # plain R integer vectors, for accuracy/AUC on the R side
  y_val_int_r   <- y_int[val_idx]

  n_train <- length(train_idx)

  best_val_loss <- Inf
  best_state <- NULL
  best_epoch <- 0
  epochs_no_improve <- 0
  history <- data.frame(epoch = integer(0),
                        train_loss = numeric(0), train_acc = numeric(0), train_auc = numeric(0),
                        valid_loss = numeric(0), valid_acc = numeric(0), valid_auc = numeric(0))

  for (ep in seq_len(epochs)) {
    net$train()
    idx <- sample.int(n_train)

    for (start in seq(1, n_train, by = batch_size)) {
      end   <- min(start + batch_size - 1, n_train)
      batch <- idx[start:end]

      opt$zero_grad()
      out  <- net(x_train[batch, , drop = FALSE])
      loss <- loss_fn(out, y_train[batch])
      loss$backward()
      opt$step()
    }

    # Report clean full-pass metrics (eval mode, dropout off) rather than
    # noisy per-batch training-mode loss -- matches the train/valid
    # loss+acc+auc reporting style used in most training scripts.
    net$eval()
    train_m <- .dl_eval_pass(net, loss_fn, x_train, y_train, y_train_int_r, n_classes)
    valid_m <- .dl_eval_pass(net, loss_fn, x_val, y_val, y_val_int_r, n_classes)

    history <- rbind(history, data.frame(
      epoch = ep,
      train_loss = train_m$loss, train_acc = train_m$acc, train_auc = train_m$auc,
      valid_loss = valid_m$loss, valid_acc = valid_m$acc, valid_auc = valid_m$auc
    ))

    if (verbose) {
      message(sprintf("Train Epoch: %d, train_loss: %.4f, train_acc: %.4f, train_auc: %.4f",
                      ep, train_m$loss, train_m$acc, train_m$auc))
      message(sprintf("Valid Epoch: %d, valid_loss: %.4f, valid_acc: %.4f, valid_auc: %.4f",
                      ep, valid_m$loss, valid_m$acc, valid_m$auc))
    }

    if (valid_m$loss < best_val_loss - 1e-6) {
      best_val_loss <- valid_m$loss
      best_epoch <- ep
      best_state <- lapply(net$state_dict(), function(t) t$clone())
      epochs_no_improve <- 0
    } else {
      epochs_no_improve <- epochs_no_improve + 1
      if (!is.null(patience) && epochs_no_improve >= patience) {
        if (verbose) message(sprintf("  early stopping at epoch %d (best epoch %d, valid_loss=%.4f)",
                                     ep, best_epoch, best_val_loss))
        break
      }
    }
  }

  if (!is.null(best_state)) {
    net$load_state_dict(best_state)
  }
  net$eval()

  list(net = net, history = history, best_epoch = best_epoch, best_val_loss = best_val_loss)
}

#' @keywords internal
#' Forward pass -> softmax probabilities as a plain R matrix (cells x classes).
.dl_predict_probs <- function(net, x_mat) {
  x_t <- torch::torch_tensor(as.matrix(x_mat), dtype = torch::torch_float())
  net$eval()
  probs <- torch::with_no_grad({
    torch::nnf_softmax(net(x_t), dim = 2)
  })
  as.matrix(torch::as_array(probs))
}


#' Train a custom deep-learning classifier with consensus label cleaning
#'
#' Drop-in deep-learning counterpart to \code{\link{train_custom_randomforest}}.
#' Same two-stage workflow: repeated stratified cross-validation is used to
#' flag unreliable labels (consensus cleaning), then a final network is
#' trained on the cleaned labelled cells. Uses a small feedforward neural
#' network (via the \code{torch} package) instead of \code{ranger}.
#'
#' @param spe A \code{SpatialExperiment} or \code{SingleCellExperiment}.
#' @param label_col Character scalar naming the label column in
#'   \code{colData(spe)} (default \code{"cutoff_label"}).
#' @param assay_name Character scalar naming the assay used as features
#'   (default \code{"norm"}).
#' @param unknown_label Character label treated as unlabeled and excluded from
#'   training/cleaning (default \code{"Unknown"}).
#' @param features Character vector of feature (marker) names to use, or
#'   \code{"all"} to use all assay rows (default \code{"all"}).
#' @param hidden_dims Integer vector of hidden layer sizes (default \code{c(64, 32)}).
#' @param dropout Dropout probability applied after each hidden layer (default \code{0.3}).
#' @param epochs Integer maximum training epochs per fold/final fit; actual
#'   training usually stops earlier via early stopping (see \code{patience})
#'   once validation loss stops improving (default \code{100}).
#' @param lr Learning rate for Adam (default \code{1e-3}).
#' @param weight_decay L2 penalty for Adam (default \code{1e-4}).
#' @param batch_size Minibatch size (default \code{256}).
#' @param cv_folds Integer number of folds used for stratified cross-validation
#'   (default \code{5}).
#' @param repeats Integer number of repeated CV rounds used during consensus
#'   cleaning (default \code{3}; kept lower than the random-forest default
#'   since each round trains a full network rather than a fast forest).
#' @param agreement_thresh Numeric in [0,1] specifying the minimum agreement rate
#'   required to keep an originally labelled cell during cleaning (default \code{0.8}).
#' @param val_frac Fraction of each training set held out (stratified by class)
#'   as a validation split for best-epoch selection and early stopping (default \code{0.15}).
#' @param patience Integer epochs to wait for validation-loss improvement before
#'   stopping early; set to \code{NULL} to disable early stopping (default \code{10}).
#' @param verbose_training Logical; if \code{TRUE}, print per-epoch train/val loss
#'   during the final model fit (the CV/cleaning-stage folds always stay quiet
#'   to avoid flooding the console). Default \code{FALSE}.
#' @param class_balanced_loss Logical; if \code{TRUE} (default), weight the loss
#'   by inverse class frequency to counter class imbalance.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A named list with components:
#' \describe{
#'   \item{spe}{The input \code{spe} with an added \code{colData} column
#'     \code{cleaned_core_label}, mirroring \code{train_custom_randomforest}.}
#'   \item{model}{An object of class \code{"cytoGateR_dl"} wrapping the trained
#'     torch network together with the feature names, class levels, and
#'     scaling stats needed to reproduce predictions.}
#'   \item{core_prob_mat}{Consensus probability matrix for the kept cells.}
#'   \item{agreement_rates}{Named numeric vector of per-cell agreement rates.}
#'   \item{features_used}{Character vector of feature names used for training.}
#'   \item{training_history}{Per-epoch train/val loss for the *final* model fit
#'     (not the CV/cleaning folds), plus which epoch was selected as best.}
#' }
#' @export
train_custom_dl <- function(spe,
                            label_col = "cutoff_label",
                            assay_name = "norm",
                            unknown_label = "Unknown",
                            features = "all",
                            hidden_dims = c(64, 32),
                            dropout = 0.3,
                            epochs = 50,
                            lr = 1e-3,
                            weight_decay = 1e-4,
                            batch_size = 256,
                            cv_folds = 5,
                            repeats = 3,
                            agreement_thresh = 0.8,
                            val_frac = 0.15,
                            patience = 10,
                            verbose_training = FALSE,
                            class_balanced_loss = TRUE,
                            seed = NULL) {



  .require_torch()
  if (!is.null(seed)) {
    set.seed(seed)
    torch::torch_manual_seed(seed)
  }
  .assert_spe(spe)



  # 1. Feature prep -- identical to train_custom_randomforest
  feat_mat_all <- SummarizedExperiment::assay(spe, assay_name)
  features_use <- if (length(features) == 1 && features == "all") {
    rownames(feat_mat_all)
  } else {
    intersect(features, rownames(feat_mat_all))
  }

  lab_vec <- SummarizedExperiment::colData(spe)[[label_col]]
  feature_df <- as.data.frame(t(feat_mat_all[features_use, , drop = FALSE]))

  core_idx <- which(lab_vec != unknown_label & !is.na(lab_vec))

  if (length(core_idx) == 0) {
    tbl <- table(lab_vec, useNA = "ifany")
    stop(
      "No core cells found: every cell in colData(spe)[['", label_col, "']] ",
      "either equals unknown_label ('", unknown_label, "') or is NA.\n",
      "Value counts for '", label_col, "':\n",
      paste(capture.output(print(tbl)), collapse = "\n"), "\n",
      "Check that 'label_col' and 'unknown_label' match how your labels are actually stored.",
      call. = FALSE
    )
  }

  core_df  <- feature_df[core_idx, , drop = FALSE]
  core_labels <- factor(lab_vec[core_idx])
  class_levels <- levels(core_labels)
  n_classes <- length(class_levels)
  rownames(core_df) <- rownames(feature_df)[core_idx]

  # 2. Scale features -- neural nets need this, ranger doesn't.
  #    Fit scaling on the full core set; reused for CV folds, final fit, and
  #    at prediction time (stored inside the returned model object).
  scale_stats <- .dl_scale_fit(as.matrix(core_df))
  core_mat_scaled <- .dl_scale_apply(as.matrix(core_df), scale_stats$center, scale_stats$scale)

  y_int <- as.integer(core_labels)  # 1-indexed: torch-for-R's cross-entropy/nll losses
  # require targets in 1..n_classes, confirmed at runtime
  # ("Indexing starts at 1 but found a 0" if 0-indexed)
  class_weights <- if (class_balanced_loss) {
    freq <- table(core_labels)
    w <- 1 / as.numeric(freq[class_levels])
    w * n_classes / sum(w)  # normalize so weights average to 1
  } else NULL

  sum_prob_mat <- matrix(0, nrow = nrow(core_df), ncol = n_classes,
                         dimnames = list(rownames(core_df), class_levels))
  match_counts <- stats::setNames(numeric(nrow(core_df)), rownames(core_df))

  # --- STAGE 1: CONSENSUS CLEANING & PROBABILITY GENERATION ---
  # Same repeated stratified CV pattern as train_custom_randomforest; the
  # ranger::ranger() fit is replaced by .dl_train_loop()/.dl_predict_probs().
  message(sprintf("Processing %d core cells via %d repeated CV folds (DL)...",
                  nrow(core_df), repeats))

  for (r in seq_len(repeats)) {
    fold_assign <- integer(nrow(core_df))
    for (cls in seq_len(n_classes)) {
      cls_idx <- which(y_int == cls)
      fold_assign[cls_idx] <- sample(rep(seq_len(cv_folds), length.out = length(cls_idx)))
    }

    for (k in seq_len(cv_folds)) {
      train_idx <- which(fold_assign != k)
      test_idx  <- which(fold_assign == k)

      fit <- .dl_train_loop(
        x_mat = core_mat_scaled[train_idx, , drop = FALSE],
        y_int = y_int[train_idx],
        n_classes = n_classes,
        hidden_dims = hidden_dims,
        dropout = dropout,
        epochs = epochs,
        lr = lr,
        weight_decay = weight_decay,
        batch_size = batch_size,
        class_weights = class_weights,
        val_frac = val_frac,
        patience = patience,
        verbose = FALSE  # CV/cleaning folds always stay quiet
      )
      net <- fit$net

      prob_preds <- .dl_predict_probs(net, core_mat_scaled[test_idx, , drop = FALSE])
      sum_prob_mat[test_idx, ] <- sum_prob_mat[test_idx, ] + prob_preds

      preds <- class_levels[max.col(prob_preds)]
      match_counts[test_idx] <- match_counts[test_idx] +
        as.numeric(preds == as.character(core_labels[test_idx]))
    }
  }

  # 3. Consensus calculations -- identical logic to the RF version
  avg_core_prob_mat <- sum_prob_mat / repeats
  agreement_rate <- match_counts / repeats

  valid_names <- names(agreement_rate)[agreement_rate >= agreement_thresh]

  new_labels <- lab_vec
  new_labels[colnames(spe) %in% names(agreement_rate[agreement_rate < agreement_thresh])] <- unknown_label
  SummarizedExperiment::colData(spe)$cleaned_core_label <- new_labels

  # 4. Final model training on cleaned cells only
  message(sprintf("Training final DL model on %d cleaned cells...", length(valid_names)))
  cleaned_mat <- core_mat_scaled[valid_names, , drop = FALSE]
  cleaned_y   <- y_int[match(valid_names, rownames(core_df))]

  final_fit <- .dl_train_loop(
    x_mat = cleaned_mat,
    y_int = cleaned_y,
    n_classes = n_classes,
    hidden_dims = hidden_dims,
    dropout = dropout,
    epochs = epochs,
    lr = lr,
    weight_decay = weight_decay,
    batch_size = batch_size,
    class_weights = class_weights,
    val_frac = val_frac,
    patience = patience,
    verbose = verbose_training
  )
  final_net <- final_fit$net
  message(sprintf("Final model: best epoch %d/%d (val_loss=%.4f)",
                  final_fit$best_epoch, nrow(final_fit$history), final_fit$best_val_loss))

  model <- structure(
    list(
      net = final_net,
      feature_names = features_use,
      class_levels = class_levels,
      center = scale_stats$center,
      scale = scale_stats$scale,
      hidden_dims = hidden_dims,
      dropout = dropout
    ),
    class = "cytoGateR_dl"
  )

  return(list(
    spe = spe,
    model = model,
    core_prob_mat = avg_core_prob_mat[valid_names, , drop = FALSE],
    agreement_rates = agreement_rate[valid_names],
    features_used = features_use,
    training_history = final_fit$history
  ))
}


#' Predict unknown cells using a trained deep-learning model
#'
#' Drop-in deep-learning counterpart to
#' \code{\link{predict_unknown_with_randomforest}}. Applies a
#' \code{"cytoGateR_dl"} model (from \code{\link{train_custom_dl}}) to cells
#' currently labelled \code{unknown_label} (or \code{NA}) and fills in
#' predictions above \code{threshold}.
#'
#' @param spe A \code{SpatialExperiment} or \code{SingleCellExperiment}.
#' @param model A \code{"cytoGateR_dl"} object, as returned by \code{train_custom_dl()$model}.
#' @param assay_name Character scalar naming the assay used as features (default \code{"norm"}).
#' @param label_col Column in \code{colData(spe)} identifying which cells are unknown (default \code{"custom_label"}).
#' @param out_col Column name to store the filled-in labels (default \code{"soft_tree_label_filled"}).
#' @param pred_col Column name to store raw predictions for the unknown subset (default \code{"dl_pred"}).
#' @param unknown_label Character label treated as unlabeled (default \code{"Unknown"}).
#' @param threshold Minimum softmax probability required to accept a prediction (default \code{0.5}).
#' @param unassigned_label Label assigned when the top probability is below \code{threshold} (default \code{"Unassigned"}).
#'
#' @return A named list with \code{spe} (updated) and \code{prob_mat} (probability
#'   matrix for the predicted-on cells only).
#' @export
predict_unknown_with_dl <- function(spe,
                                    model,
                                    assay_name = "norm",
                                    label_col = "custom_label",
                                    out_col = "soft_tree_label_filled",
                                    pred_col = "dl_pred",
                                    unknown_label = "Unknown",
                                    threshold = 0.5,
                                    unassigned_label = "Unassigned") {
  .require_torch()
  .assert_spe(spe)
  if (!inherits(model, "cytoGateR_dl")) stop("model must be a 'cytoGateR_dl' object.")

  # 1. Identify target cells
  labels <- SummarizedExperiment::colData(spe)[[label_col]]
  replace_idx <- is.na(labels) | labels == unknown_label

  if (sum(replace_idx) == 0) {
    message("No unknown cells found to predict.")
    return(list(spe = spe, prob_mat = NULL))
  }

  # 2. Prepare features ONLY for unknown cells, aligned to training feature order
  feat_mat <- SummarizedExperiment::assay(spe, assay_name)
  feature_df <- as.data.frame(t(feat_mat[, replace_idx, drop = FALSE]))
  feature_df <- feature_df[, model$feature_names, drop = FALSE]

  scaled_mat <- .dl_scale_apply(as.matrix(feature_df), model$center, model$scale)

  # 3. Predict ONLY on unknowns
  prob_mat <- .dl_predict_probs(model$net, scaled_mat)
  dimnames(prob_mat) <- list(rownames(feature_df), model$class_levels)

  # 4. Determine winners and confidence for the subset
  max_probs <- apply(prob_mat, 1, max)
  winning_indices <- apply(prob_mat, 1, which.max)
  raw_preds <- colnames(prob_mat)[winning_indices]

  final_preds_subset <- ifelse(max_probs >= threshold, raw_preds, unassigned_label)

  # 5. Fill back into the full SPE -- identical bookkeeping to the RF version
  full_final_preds <- rep(NA_character_, ncol(spe))
  full_confidence <- rep(NA_character_, ncol(spe))
  filled <- labels

  full_final_preds[replace_idx] <- final_preds_subset
  full_confidence[replace_idx] <- max_probs
  filled[replace_idx] <- final_preds_subset

  SummarizedExperiment::colData(spe)[[pred_col]] <- full_final_preds
  SummarizedExperiment::colData(spe)[[out_col]] <- filled
  SummarizedExperiment::colData(spe)$dl_confidence <- as.numeric(full_confidence)

  return(list(spe = spe, prob_mat = prob_mat))
}
