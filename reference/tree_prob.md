# Tree probability for one cell

Compute the probability that a single cell belongs to a given cell type
by combining per-node scores along a tree path, applying a
complexity-based weight, and an optional negative marker penalty.

## Usage

``` r
tree_prob(
  tree,
  expr_mat,
  cell_i,
  combine = c("mean", "product"),
  neg_markers = character(0),
  marker_stats = NULL,
  neg_strength = 0.8,
  lambda = 2
)
```

## Arguments

- tree:

  A gating tree as produced by \[build_fullcoverage_tree()\].

- expr_mat:

  Numeric matrix-like expression object with markers in rows and cells
  in columns.

- cell_i:

  Integer index of the cell/column to score.

- combine:

  How to combine node-wise scores along the path. One of \`"mean"\` or
  \`"product"\`.

- neg_markers:

  Character vector of negative marker names used for penalty.

- marker_stats:

  A list of marker statistics (including cutoff and scale) produced by
  \[fit_marker_stats()\] used to calculate the penalty.

- neg_strength:

  Numeric scalar in \[0,1\] defining how strongly negative marker
  expression should penalize the final score. Default is 0.8.

- lambda:

  Numeric scalar representing the "Skepticism Factor." Higher values
  place a heavier penalty on cell types defined by fewer markers.
  Default is 0.2.

## Value

Numeric scalar probability (typically in \[0,1\]) or \`NA_real\_\` if no
usable node scores are available for the cell.

## Examples

``` r
tree <- list(
  type = "node", marker = "CD3", cutoff = 1, scale = 0.5,
  left = list(type = "leaf"), right = list(type = "leaf")
)
expr_mat <- rbind(CD3 = c(0.5, 1.5))
tree_prob(tree, expr_mat, cell_i = 2)
#> [1] 0.2436862
```
