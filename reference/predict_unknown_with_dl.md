# Predict unknown cells using a trained deep-learning model

Drop-in deep-learning counterpart to
[`predict_unknown_with_randomforest`](https://sydneybiox.github.io/cytoGateR/reference/predict_unknown_with_randomforest.md).
Applies a `"cytoGateR_dl"` model (from
[`train_custom_dl`](https://sydneybiox.github.io/cytoGateR/reference/train_custom_dl.md))
to cells currently labelled `unknown_label` (or `NA`) and fills in
predictions above `threshold`.

## Usage

``` r
predict_unknown_with_dl(
  spe,
  model,
  assay_name = "exprs",
  label_col = "custom_label",
  out_col = "soft_tree_label_filled",
  pred_col = "dl_pred",
  unknown_label = "Unknown",
  threshold = 0.5,
  unassigned_label = "Unassigned"
)
```

## Arguments

- spe:

  A `SpatialExperiment` or `SingleCellExperiment`.

- model:

  A `"cytoGateR_dl"` object, as returned by `train_custom_dl()$model`.

- assay_name:

  Character scalar naming the assay used as features (default
  `"exprs"`).

- label_col:

  Column in `colData(spe)` identifying which cells are unknown (default
  `"custom_label"`).

- out_col:

  Column name to store the filled-in labels (default
  `"soft_tree_label_filled"`).

- pred_col:

  Column name to store raw predictions for the unknown subset (default
  `"dl_pred"`).

- unknown_label:

  Character label treated as unlabeled (default `"Unknown"`).

- threshold:

  Minimum softmax probability required to accept a prediction (default
  `0.5`).

- unassigned_label:

  Label assigned when the top probability is below `threshold` (default
  `"Unassigned"`).

## Value

A named list with `spe` (updated) and `prob_mat` (probability matrix for
the predicted-on cells only).

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
  dl_prediction <- predict_unknown_with_dl(
    dl_fit$spe,
    dl_fit$model,
    label_col = "cleaned_core_label",
    unknown_label = "Unassigned",
    threshold = 0.5
  )
  head(SummarizedExperiment::colData(dl_prediction$spe)$dl_pred)
}
# }
```
