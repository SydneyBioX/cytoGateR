# Run soft-gating pipeline end-to-end

Run soft-gating pipeline end-to-end

## Usage

``` r
run_soft_gating(
  spe,
  lineage_table,
  assay_name = "exprs",
  unknown_thresh = 0.4,
  store = TRUE
)
```

## Arguments

- spe:

  SpatialExperiment / SingleCellExperiment.

- lineage_table:

  Tibble with cell_type/pos_markers/neg_markers.

- assay_name:

  Assay to use (default "exprs").

- unknown_thresh:

  Threshold for Unknown label (default 0.4).

- store:

  Store results back into spe (default TRUE).

## Value

A list with: spe (if store=TRUE), marker_stats, prob_mat, labels.

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
soft_result <- run_soft_gating(spe, lineage_table)
head(soft_result$labels)
#> [1] "Tcell"   "Bcell"   "Unknown" "Unknown" "Unknown" "Tcell"  
```
