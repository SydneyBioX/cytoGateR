# Build full-coverage gating tree

Build full-coverage gating tree

## Usage

``` r
build_fullcoverage_tree(
  expr_mat,
  markers_pos,
  markers_neg = character(0),
  cell_idx = seq_len(ncol(expr_mat)),
  depth = 0,
  max_depth = 4,
  min_cells = 200,
  min_score = 0.5,
  cutoff_method = c("mean", "equal_posteriors"),
  gmm_model_names = NULL
)
```

## Arguments

- expr_mat:

  expression matrix

- markers_pos:

  positive markers

- markers_neg:

  negative markers

- cell_idx:

  indices

- depth:

  depth

- max_depth:

  max depth

- min_cells:

  minimum cells

- min_score:

  minimum separability

- cutoff_method:

  Cutoff selection method passed to \[fit_gmm_2()\].

- gmm_model_names:

  Optional character vector of model names to pass to
  \[mclust::Mclust()\] (e.g., "V" to forbid equal-variance in 1D).

## Value

tree object

## Examples

``` r
set.seed(1)
expr_mat <- rbind(
  CD3 = c(rnorm(50, 0, 0.2), rnorm(50, 2, 0.2)),
  CD20 = c(rnorm(50, 2, 0.2), rnorm(50, 0, 0.2))
)
tree <- build_fullcoverage_tree(
  expr_mat,
  markers_pos = c("CD3", "CD20"),
  max_depth = 2,
  min_cells = 20,
  min_score = 0.1
)
tree$type
#> [1] "node"
```
