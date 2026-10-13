# Calculate F1 Score per Cell Type

Calculate F1 Score per Cell Type

## Usage

``` r
calculate_f1(spe, ref_col = "ref_broad", pred_col = "pred_broad")
```

## Arguments

- spe:

  SpatialExperiment object

- ref_col:

  Column name for ground truth (e.g., "manual_label")

- pred_col:

  Column name for results (e.g., "hier_label")

## Value

A data frame with Precision, Recall, and F1 score per type

## Examples

``` r
data("cytoGateR_example", package = "cytoGateR")
image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
spe <- cytoGateR_example[, image_ids == image_ids[1L]]
calculate_f1(spe, ref_col = "cell_type", pred_col = "knn_label")
#>       Category Precision Recall F1_Score Cell_Count
#> 1          CD8     0.774  0.568    0.655        736
#> 2        Bcell     0.450  0.880    0.596        284
#> 3      BnTcell     0.654  0.830    0.732       1033
#> 4          CD4     0.724  0.751    0.737        747
#> 5         Treg     0.446  0.565    0.498        138
#> 6      Myeloid     0.899  0.403    0.556        683
#> 7       Stroma     0.770  0.694    0.730        444
#> 8   Neutrophil     0.341  0.667    0.452         21
#> 9  Plasma_cell     0.246  0.556    0.341         27
#> 10       Tumor     0.906  0.875    0.890         88
```
