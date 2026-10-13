# Add custom labels to a SpatialExperiment

Applies the confidence-based labelling scheme to a probability matrix
and stores the resulting labels in \`colData(spe)\`.

## Usage

``` r
custom_labels(
  spe,
  prob_mat,
  colname = "custom_label",
  base_thresh = 0.4,
  quantile_cut = 0.75,
  iqr_mult = 1.5,
  unknown_label = "Unknown",
  flag_fn = NULL
)
```

## Arguments

- spe:

  A \`SpatialExperiment\` or \`SingleCellExperiment\`.

- prob_mat:

  Numeric probability matrix (cells x lineages). Rows must correspond to
  \`colnames(spe)\`; if rownames are present they are checked.

- colname:

  Name of the output column in \`colData\` (default "custom_label").

- base_thresh:

  Minimum allowable threshold (default 0.4).

- quantile_cut:

  Quantile used to set the adaptive threshold (default 0.75).

- iqr_mult:

  Multiplier applied to the IQR (default 1.5).

- unknown_label:

  Label assigned to uncertain cells (default "Unknown").

- flag_fn:

  Optional function taking a probability vector and returning a logical
  indicating confidence; defaults to the internal Tukey-style rule.

## Value

The input \`spe\` with a new \`colData\` column containing custom
labels.

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
spe <- custom_labels(
  spe,
  prob_mat,
  flag_fn = function(x) max(x) >= 0.8
)
table(SummarizedExperiment::colData(spe)$custom_label)
#> 
#>       Bcell     BnTcell         CD4         CD8     Myeloid  Neutrophil 
#>         401         976         505         381         210          30 
#> Plasma_cell      Stroma        Treg       Tumor     Unknown 
#>          17         347         110          70        1316 
```
