# Plot a marker priority "tree" for a given cell type (weights order)

This is a simple visualization: markers ordered by weight and connected
as a chain ROOT -\> marker1 -\> marker2 -\> ...

## Usage

``` r
plot_marker_priority_tree(lineage_table, marker_stats, cell_type_name)
```

## Arguments

- lineage_table:

  Tibble with cell_type/pos_markers/neg_markers.

- marker_stats:

  Output of fit_marker_stats().

- cell_type_name:

  Which cell type to plot.

## Value

A ggraph/ggplot object (requires igraph + ggraph).

## Examples

``` r
markers <- c("CD3e", "CD20")
lineage <- data.frame(
  cell_type = "Immune",
  pos_markers = I(list(markers)),
  neg_markers = I(list(character()))
)
marker_stats <- list(
  CD3e = list(weight = 0.9),
  CD20 = list(weight = 0.7)
)
plot_marker_priority_tree(lineage, marker_stats, "Immune")
```
