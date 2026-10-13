# Apply Protected Spatial Priors to kNN Probabilities

Apply Protected Spatial Priors to kNN Probabilities

## Usage

``` r
calculate_spatial_prior_labels(
  spe,
  prob_mat,
  k_spatial = 50,
  lambda = 0.2,
  protect_threshold = 0.85,
  out_col = "knn_spatial_label"
)
```

## Arguments

- spe:

  SpatialExperiment object.

- prob_mat:

  The full probability matrix \[Cells x Types\] from
  predict_unknown_with_knn.

- k_spatial:

  Number of physical neighbors for the "Prior" (default 50).

- lambda:

  The "Spatial Weight" (default 0.2).

- protect_threshold:

  Confidence level above which spatial priors are ignored (default
  0.85).

- out_col:

  Column name for the new labels.

## Value

The input \`SpatialExperiment\` with the spatially adjusted labels
stored in \`colData()\` under \`out_col\`.

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
spe <- calculate_spatial_prior_labels(
  spe,
  prob_mat,
  k_spatial = 3,
  out_col = "spatial_label"
)
#> Integrating spatial context (Protecting cells > 0.85 confidence)...
table(SummarizedExperiment::colData(spe)$spatial_label)
#> 
#>       Bcell     BnTcell         CD4         CD8     Myeloid  Neutrophil 
#>         714        1392         777         514         274          36 
#> Plasma_cell      Stroma        Treg       Tumor 
#>          40         398         135          83 
```
