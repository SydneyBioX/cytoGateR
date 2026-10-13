# Fit per-marker 2-GMM statistics for soft gating

Fits a 2-component Gaussian Mixture Model (GMM) per marker and returns a
list of per-marker stats (cutoff, scale, weight, mu_high, valid).

## Usage

``` r
fit_marker_stats(spe, markers, assay_name = "exprs", min_n = 50L)
```

## Arguments

- spe:

  SpatialExperiment / SingleCellExperiment.

- markers:

  Character vector of marker names.

- assay_name:

  Assay to use (default "exprs").

- min_n:

  Minimum finite observations required to fit GMM (default 50).

## Value

Named list. Each element is a list with fields: valid, cutoff, scale,
weight, mu_high.

## Examples

``` r
data("cytoGateR_example", package = "cytoGateR")
image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
spe <- cytoGateR_example[, image_ids == image_ids[1L]]
marker_stats <- fit_marker_stats(
  spe,
  markers = c("CD3e", "CD20")
)
marker_stats
#> $CD3e
#> $CD3e$valid
#> [1] TRUE
#> 
#> $CD3e$cutoff
#> [1] 1.565636
#> 
#> $CD3e$scale
#> [1] 0.4582558
#> 
#> $CD3e$weight
#> [1] 2.783092
#> 
#> $CD3e$mu_high
#> [1] 2.482352
#> 
#> 
#> $CD20
#> $CD20$valid
#> [1] TRUE
#> 
#> $CD20$cutoff
#> [1] 1.188319
#> 
#> $CD20$scale
#> [1] 0.5642736
#> 
#> $CD20$weight
#> [1] 1.718334
#> 
#> $CD20$mu_high
#> [1] 2.015195
#> 
#> 
```
