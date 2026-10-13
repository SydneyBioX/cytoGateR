# Build Lineage Hierarchy from Core Cells

Build Lineage Hierarchy from Core Cells

## Usage

``` r
build_lineage_hierarchy(
  spe,
  label_col = "cleaned_core_label",
  assay_name = "exprs"
)
```

## Arguments

- spe:

  A SpatialExperiment object with gated core labels.

- label_col:

  The column containing gated labels (e.g., "cleaned_core_label").

- assay_name:

  Assay to use for pseudobulk calculation.

## Value

An hclust object representing the cell type hierarchy.

## Examples

``` r
data("cytoGateR_example", package = "cytoGateR")
image_ids <- SummarizedExperiment::colData(cytoGateR_example)$image_name
spe <- cytoGateR_example[, image_ids == image_ids[1L]]
SummarizedExperiment::colData(spe)$cleaned_core_label <-
  SummarizedExperiment::colData(spe)$core_group
SummarizedExperiment::colData(spe)$cleaned_core_label[
  SummarizedExperiment::colData(spe)$cleaned_core_label == "Unassigned"
] <- "Unknown"
hc_tree <- build_lineage_hierarchy(spe)
hc_tree$labels
#>  [1] "BnTcell"     "CD4"         "Bcell"       "Treg"        "Plasma_cell"
#>  [6] "Myeloid"     "CD8"         "Stroma"      "Tumor"       "Neutrophil" 
plot(hc_tree)
```
