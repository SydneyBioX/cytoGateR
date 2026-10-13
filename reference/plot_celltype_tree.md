# Plot a cell-type gating tree

Plot a cell-type gating tree

## Usage

``` r
plot_celltype_tree(tree, title = "")
```

## Arguments

- tree:

  A gating tree produced by \[build_fullcoverage_tree()\].

- title:

  Plot title string.

## Value

A \`ggplot\` object.

## Examples

``` r
tree <- list(
  type = "node", depth = 0, cells = 1:6, marker = "CD3e",
  cutoff = 0.5, sep_score = 0.8,
  left = list(type = "leaf", depth = 1, cells = 1:3),
  right = list(type = "leaf", depth = 1, cells = 4:6)
)
plot_celltype_tree(tree, title = "Example gating tree")
```
