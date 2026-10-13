# Apply cutoff-based labels to a result list

Takes a \`res\` object (as returned by \`run_soft_gating()\` or
\`run_tree_gating()\`), builds a per-cell-type label matrix using a
cutoff function, stores it as \`res\$label_mat\`, and writes a
single-label assignment into \`colData(res\$spe)\`.

## Usage

``` r
apply_cutoff_labels(
  res,
  cutoff_fn = prob_quantile_cutoff,
  label_col = "cutoff_label",
  unknown_label = "Unknown"
)
```

## Arguments

- res:

  List containing at least \`spe\` and \`prob_mat\`.

- cutoff_fn:

  Function passed to \[probability_label_matrix()\] to compute
  per-column cutoffs (default \[prob_quantile_cutoff()\]).

- label_col:

  Name of the column to add to \`colData(res\$spe)\` (default
  "cutoff_label").

- unknown_label:

  Label used when zero or multiple cell types pass their cutoff (default
  "Unknown").

## Value

The input \`res\` list with updated \`spe\` and a new \`label_mat\`
element.

## Details

A cell receives a label only if exactly one cell type passes its cutoff;
if multiple (or zero) cell types pass, the cell is labeled
\`unknown_label\`.

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
result <- apply_cutoff_labels(
  list(spe = spe, prob_mat = prob_mat),
  cutoff_fn = function(x) 0.5
)
table(SummarizedExperiment::colData(result$spe)$cutoff_label)
#> 
#>       Bcell     BnTcell         CD4         CD8     Myeloid  Neutrophil 
#>         554        1310         775         540         306          41 
#> Plasma_cell      Stroma        Treg       Tumor     Unknown 
#>          60         400         175          85         117 
```
