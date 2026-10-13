# Build Hierarchical Reference with Flexible Marker Selection

Build Hierarchical Reference with Flexible Marker Selection

## Usage

``` r
build_hierarchical_reference(
  spe,
  hc_tree,
  marker_stats = NULL,
  label_col = "cleaned_core_label",
  top_n = 5,
  assay_name = "exprs",
  unknown_label = "Unassigned"
)
```

## Arguments

- spe:

  SpatialExperiment containing core cells.

- hc_tree:

  The hclust object from build_lineage_hierarchy.

- marker_stats:

  Output from fit_marker_stats(). Required if top_n is not NULL.

- label_col:

  Column with core labels (e.g., "cleaned_core_label").

- top_n:

  Number of top DE markers to use per node. Set to NULL to use all
  markers.

- assay_name:

  Assay to use (default "exprs").

- unknown_label:

  Character label for unassigned cells (default "Unassigned").

## Value

A named list of node-specific references. Each element contains training
data, binary left/right labels, selected markers, and the lineage
members represented by each side of the node.

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
hier_ref <- build_hierarchical_reference(
  spe,
  hc_tree,
  top_n = NULL,
  unknown_label = "Unknown"
)
names(hier_ref)
#> [1] "Node_1" "Node_2" "Node_3" "Node_4" "Node_5" "Node_6" "Node_7" "Node_8"
#> [9] "Node_9"
hier_ref[[1L]]$markers
#>  [1] "MPO"    "H3"     "SMA"    "CD16"   "CD38"   "HLA_DR" "CD27"   "CD15"  
#>  [9] "CD45RA" "CD163"  "B2M"    "CD20"   "CD68"   "IDO1"   "CD3e"   "LAG3"  
#> [17] "CD11c"  "PD_1"   "PDGFRB" "CD7"    "GZMB"   "PD_L1"  "TCF7"   "CD45RO"
#> [25] "FOXP3"  "ICOS"   "CD8a"   "CA9"    "CD33"   "Ki67"   "VISTA"  "CD40"  
#> [33] "CD4"    "CD14"   "CDH1"   "CD303"  "CD206"  "c_PARP" "DNA1"   "DNA2"  
```
