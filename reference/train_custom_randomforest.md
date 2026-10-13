# Train a custom random forest with consensus label cleaning

Train a `ranger` random forest classifier on labelled cells from a
`SpatialExperiment` or `SingleCellExperiment`. The workflow has two
stages:

1.  **Consensus cleaning:** repeated stratified cross-validation is used
    to estimate how consistently each labelled cell is predicted as its
    original label. Cells with agreement below `agreement_thresh` are
    flagged as unreliable and relabeled as `unknown_label` in
    `colData(spe)$cleaned_core_label`.

2.  **Model evaluation and final training:** cross-validation is run on
    the cleaned labelled set to report accuracy and macro-F1, and a
    final probability random forest model is trained on all cleaned
    labelled cells.

## Usage

``` r
train_custom_randomforest(
  spe,
  label_col = "cutoff_label",
  assay_name = "exprs",
  unknown_label = "Unknown",
  num.trees = 200,
  mtry = NULL,
  seed = NULL,
  features = "all",
  cv_folds = 5,
  repeats = 10,
  agreement_thresh = 0.8,
  num_threads = 2
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

- num.trees:

  Integer number of trees for the final `ranger` model (default `200`).

- mtry:

  Optional integer `mtry` for `ranger`. If `NULL`, uses `floor(sqrt(p))`
  where `p` is the number of features used.

- seed:

  Optional integer seed for reproducibility.

- features:

  Character vector of feature (marker) names to use, or `"all"` to use
  all assay rows (default `"all"`).

- cv_folds:

  Integer number of folds used for stratified cross-validation (default
  `5`).

- repeats:

  Integer number of repeated CV rounds used during consensus cleaning
  (default `10`).

- agreement_thresh:

  Numeric in \[0,1\] specifying the minimum agreement rate required to
  keep an originally labelled cell during cleaning (default `0.8`).

- num_threads:

  Integer requested number of threads for model fitting.

## Value

A named list with components:

- spe:

  The input `spe` with an added `colData` column `cleaned_core_label`,
  where low-agreement labelled cells are set to `unknown_label`.

- model:

  A fitted `ranger` model trained on the cleaned labelled cells with
  `probability = TRUE`.

- agreement_rates:

  Named numeric vector of per-cell agreement rates from the consensus
  cleaning stage (indexed by training-row names).

- features_used:

  Character vector of feature names used for training.

- metrics:

  List of evaluation outputs on the cleaned labelled set:

  confusion_matrix

  :   Confusion matrix for pooled CV predictions.

  accuracy

  :   Overall pooled CV accuracy.

  f1_macro

  :   Overall pooled macro-F1.

  cv_overall

  :   Per-fold data.frame with accuracy and macro-F1.

  cv_class_metrics

  :   Per-fold, per-class precision/recall/F1 table.

- test_pred:

  Factor of pooled cross-validated predictions on cleaned labelled
  cells.

- test_truth:

  Factor of pooled cross-validated true labels on cleaned labelled
  cells.

## Details

During consensus cleaning, models are trained without probability
estimation. The final model is trained with `probability = TRUE` to
support downstream thresholding workflows.

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
fit$model
#> Ranger result
#> 
#> Call:
#>  ranger::ranger(original_label ~ ., data = cleaned_df, num.trees = num.trees,      mtry = mtry_val, probability = TRUE, num.threads = num_threads) 
#> 
#> Type:                             Probability estimation 
#> Number of trees:                  25 
#> Sample size:                      1378 
#> Number of independent variables:  4 
#> Mtry:                             2 
#> Target node size:                 10 
#> Variable importance mode:         none 
#> Splitrule:                        gini 
#> OOB prediction error (Brier s.):  0.1906214 
# }
```
