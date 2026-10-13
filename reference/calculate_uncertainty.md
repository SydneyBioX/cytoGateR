# Calculate Full Uncertainty Suite with Sample Grouping

Calculate Full Uncertainty Suite with Sample Grouping

## Usage

``` r
calculate_uncertainty(prob_mat, spe, sample_col = "sample_id", k_spatial = 15)
```

## Arguments

- prob_mat:

  The \[Cells x Types\] matrix from your kNN or randomforest

- spe:

  The SpatialExperiment object.

- sample_col:

  The column identifying different images/samples.

- k_spatial:

  Neighbors for the spatial discordance check (default 15). from
  training) to rbind with prob_mat before processing. Both must share
  the same column names.

## Value

A data frame with one row per cell and columns for normalized entropy,
Gini impurity, margin uncertainty, and spatial discordance.

## Examples

``` r
data("cytoGateR_example", package = "cytoGateR")
image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
spe <- cytoGateR_example[, image_ids == image_ids[1L]]
prob_cols <- grep(
  "^KNN_P_",
  names(SummarizedExperiment::colData(spe)),
  value = TRUE
)
prob_mat <- as.matrix(SummarizedExperiment::colData(spe)[, prob_cols])
rownames(prob_mat) <- colnames(spe)
colnames(prob_mat) <- sub("^KNN_P_", "", colnames(prob_mat))
uncertainty <- calculate_uncertainty(
  prob_mat,
  spe,
  sample_col = "sample_id",
  k_spatial = 3
)
head(uncertainty)
#>      cell_id entropy gini_impurity margin_uncertainty spatial_discordance
#> 12_1    12_1       0             0                  0           1.0000000
#> 12_2    12_2       0             0                  0           0.3333333
#> 12_3    12_3       0             0                  0           0.0000000
#> 12_4    12_4       0             0                  0           0.6666667
#> 12_5    12_5       0             0                  0           0.6666667
#> 12_6    12_6       0             0                  0           0.6666667
```
