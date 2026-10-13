# Calculate soft gating probabilities per cell type

Computes a probability matrix \[cells x cell_types\] by combining
per-marker logistic scores weighted by separability, with optional
negative-marker penalties.

## Usage

``` r
calculate_soft_scores(
  spe,
  marker_stats,
  lineage_table,
  assay_name = "exprs",
  neg_strength = 0.8
)
```

## Arguments

- spe:

  SpatialExperiment / SingleCellExperiment.

- marker_stats:

  Output of fit_marker_stats().

- lineage_table:

  Tibble with columns: cell_type, pos_markers (list), neg_markers
  (list).

- assay_name:

  Assay to use (default "exprs").

- neg_strength:

  Penalty strength multiplier for negative markers (default 0.8).

## Value

Numeric matrix \[ncol(spe) x n_types\].

## Examples

``` r
data("cytoGateR_example", package = "cytoGateR")
image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
spe <- cytoGateR_example[, image_ids == image_ids[1L]]
lineage_table <- data.frame(
  cell_type = c("Tcell", "Bcell"),
  pos_markers = I(list("CD3e", "CD20")),
  neg_markers = I(list("CD20", "CD3e"))
)
marker_stats <- fit_marker_stats(
  spe,
  markers = c("CD3e", "CD20")
)
prob_mat <- calculate_soft_scores(spe, marker_stats, lineage_table)
head(prob_mat)
#>          Tcell      Bcell
#> 12_1 0.5396101 0.12148770
#> 12_2 0.1372688 0.54922936
#> 12_3 0.2033217 0.21258899
#> 12_4 0.3451912 0.16985400
#> 12_5 0.3092308 0.17932174
#> 12_6 0.7137421 0.07463277
```
