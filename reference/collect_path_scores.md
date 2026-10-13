# Collect scores along a tree path

Traverse a gating tree for a single cell and compute node-wise marker
scores along the deterministic path defined by comparing marker
expression to each node's cutoff.

## Usage

``` r
collect_path_scores(tree, expr_mat, cell_i)
```

## Arguments

- tree:

  A gating tree as produced by \[build_fullcoverage_tree()\].

- expr_mat:

  Numeric matrix-like expression object with markers in rows and cells
  in columns (e.g. \`assay(spe, "norm")\`). Row names must include
  marker names.

- cell_i:

  Integer index of the cell/column to score.

## Value

A numeric vector of per-node scores (may contain \`NA_real\_\` if a
marker value is not finite for that cell).

## Examples

``` r
tree <- list(
  type = "node", marker = "CD3", cutoff = 1, scale = 0.5,
  left = list(type = "leaf"), right = list(type = "leaf")
)
expr_mat <- rbind(CD3 = c(0.5, 1.5))
collect_path_scores(tree, expr_mat, cell_i = 1)
#> [1] 0.2689414
```
