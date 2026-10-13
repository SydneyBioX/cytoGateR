# Train a custom deep-learning classifier with consensus label cleaning

Drop-in deep-learning counterpart to
[`train_custom_randomforest`](https://sydneybiox.github.io/cytoGateR/reference/train_custom_randomforest.md).
Same two-stage workflow: repeated stratified cross-validation is used to
flag unreliable labels (consensus cleaning), then a final network is
trained on the cleaned labelled cells. Uses a small feedforward neural
network (via the `torch` package) instead of `ranger`.

## Usage

``` r
train_custom_dl(
  spe,
  label_col = "cutoff_label",
  assay_name = "exprs",
  unknown_label = "Unknown",
  features = "all",
  hidden_dims = c(64, 32),
  dropout = 0.3,
  epochs = 50,
  lr = 0.001,
  weight_decay = 1e-04,
  batch_size = 256,
  cv_folds = 5,
  repeats = 3,
  agreement_thresh = 0.8,
  val_frac = 0.15,
  patience = 10,
  verbose_training = FALSE,
  class_balanced_loss = TRUE,
  seed = NULL
)
```

## Arguments

- spe:

  A `SpatialExperiment` or `SingleCellExperiment`.

- label_col:

  Character scalar naming the label column in `colData(spe)` (default
  `"cutoff_label"`).

- assay_name:

  Character scalar naming the assay used as features (default
  `"exprs"`).

- unknown_label:

  Character label treated as unlabeled and excluded from
  training/cleaning (default `"Unknown"`).

- features:

  Character vector of feature (marker) names to use, or `"all"` to use
  all assay rows (default `"all"`).

- hidden_dims:

  Integer vector of hidden layer sizes (default `c(64, 32)`).

- dropout:

  Dropout probability applied after each hidden layer (default `0.3`).

- epochs:

  Integer maximum training epochs per fold/final fit; actual training
  usually stops earlier via early stopping (see `patience`) once
  validation loss stops improving (default `100`).

- lr:

  Learning rate for Adam (default `1e-3`).

- weight_decay:

  L2 penalty for Adam (default `1e-4`).

- batch_size:

  Minibatch size (default `256`).

- cv_folds:

  Integer number of folds used for stratified cross-validation (default
  `5`).

- repeats:

  Integer number of repeated CV rounds used during consensus cleaning
  (default `3`; kept lower than the random-forest default since each
  round trains a full network rather than a fast forest).

- agreement_thresh:

  Numeric in \[0,1\] specifying the minimum agreement rate required to
  keep an originally labelled cell during cleaning (default `0.8`).

- val_frac:

  Fraction of each training set held out (stratified by class) as a
  validation split for best-epoch selection and early stopping (default
  `0.15`).

- patience:

  Integer epochs to wait for validation-loss improvement before stopping
  early; set to `NULL` to disable early stopping (default `10`).

- verbose_training:

  Logical; if `TRUE`, print per-epoch train/val loss during the final
  model fit (the CV/cleaning-stage folds always stay quiet to avoid
  flooding the console). Default `FALSE`.

- class_balanced_loss:

  Logical; if `TRUE` (default), weight the loss by inverse class
  frequency to counter class imbalance.

- seed:

  Optional integer seed for reproducibility.

## Value

A named list with components:

- spe:

  The input `spe` with an added `colData` column `cleaned_core_label`,
  mirroring `train_custom_randomforest`.

- model:

  An object of class `"cytoGateR_dl"` wrapping the trained torch network
  together with the feature names, class levels, and scaling stats
  needed to reproduce predictions.

- core_prob_mat:

  Consensus probability matrix for the kept cells.

- agreement_rates:

  Named numeric vector of per-cell agreement rates.

- features_used:

  Character vector of feature names used for training.

- training_history:

  Per-epoch train/val loss for the \*final\* model fit (not the
  CV/cleaning folds), plus which epoch was selected as best.

## Examples

``` r
# \donttest{
if (requireNamespace("torch", quietly = TRUE) &&
    torch::torch_is_installed()) {
  data("cytoGateR_example", package = "cytoGateR")
  image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
  spe <- cytoGateR_example[, image_ids == image_ids[1L]]
  dl_fit <- train_custom_dl(
    spe,
    label_col = "core_group",
    unknown_label = "Unassigned",
    features = c("CD3e", "CD20", "CD4", "CD8a"),
    hidden_dims = 8,
    epochs = 2,
    batch_size = 256,
    cv_folds = 2,
    repeats = 1,
    agreement_thresh = 0,
    patience = 1,
    seed = 1
  )
  dl_fit$training_history
}
# }
```
