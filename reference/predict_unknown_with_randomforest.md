# Predict and fill unknown labels using a trained random forest model

Use a trained `ranger` classification model to predict cell-type labels
for all cells in a `SpatialExperiment` or `SingleCellExperiment` object.
Predicted labels are used to fill entries in `label_col` that are either
`NA` or equal to `unknown_label`.

## Usage

``` r
predict_unknown_with_randomforest(
  spe,
  model,
  assay_name = "exprs",
  label_col = "custom_label",
  out_col = "soft_tree_label_filled",
  pred_col = "rf_pred",
  unknown_label = "Unknown",
  threshold = 0.5,
  unassigned_label = "Unassigned"
)
```

## Arguments

- spe:

  A `SpatialExperiment` or `SingleCellExperiment` object.

- model:

  A fitted `ranger` model returned by `train_custom_rf()`. The model
  should be trained with `probability = TRUE` to enable threshold-based
  filtering.

- assay_name:

  Character scalar specifying which assay to use as feature input
  (default `"exprs"`).

- label_col:

  Character scalar naming the existing label column in `colData(spe)`
  (default `"custom_label"`).

- out_col:

  Character scalar naming the output column that will contain the filled
  labels (default `"soft_tree_label_filled"`).

- pred_col:

  Character scalar naming the column used to store raw model predictions
  (default `"rf_pred"`).

- unknown_label:

  Character value treated as missing and eligible for replacement
  (default `"Unknown"`).

- threshold:

  Numeric value in \[0,1\] specifying the minimum required class
  probability for a prediction to be accepted (default `0.5`). Ignored
  if the model was not trained with `probability = TRUE`.

- unassigned_label:

  Character label assigned to cells whose maximum class probability is
  below `threshold` (default `"Unassigned"`).

## Value

The input `spe` object with updated `colData` columns:

- `pred_col`:

  Raw predicted labels from the random forest.

- `out_col`:

  Filled labels, replacing only `NA` or `unknown_label` entries in
  `label_col`.

- `rf_confidence`:

  Maximum predicted class probability per cell.

## Details

If the model was trained with `probability = TRUE`, class probabilities
are used to compute a confidence score for each prediction. Predictions
with maximum class probability below `threshold` are reassigned to
`unassigned_label`. The maximum probability is stored for quality
control.

Feature columns are automatically aligned to match the variables used
during model training. Only cells originally labeled as `NA` or
`unknown_label` are replaced in the output column.

## Examples

``` r
# \donttest{
data("cytoGateR_example", package = "cytoGateR")
image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
spe <- cytoGateR_example[, image_ids == image_ids[1L]]
fit <- train_custom_randomforest(
  spe,
  label_col = "core_group",
  unknown_label = "Unassigned",
  features = c("CD3e", "CD20", "CD4", "CD8a"),
  num.trees = 25,
  cv_folds = 2,
  repeats = 1,
  agreement_thresh = 0,
  num_threads = 1,
  seed = 1
)
#> Processing 1378 core cells via 1 repeated CV folds...
#> Training final model on 1378 cleaned cells...
prediction <- predict_unknown_with_randomforest(
  fit$spe,
  fit$model,
  label_col = "cleaned_core_label",
  unknown_label = "Unassigned",
  threshold = 0.5
)
table(SummarizedExperiment::colData(prediction$spe)$soft_tree_label_filled)
#> 
#>       Bcell     BnTcell         CD4         CD8     Myeloid  Neutrophil 
#>         578        1389         422         391         230          12 
#> Plasma_cell      Stroma        Treg       Tumor  Unassigned 
#>          41         330         230          27         713 
# }
```
