# Build a binary label matrix from probabilities

Applies a per-cell-type cutoff function to each probability column to
build a logical matrix indicating which cell types pass their cutoff for
each cell.

## Usage

``` r
probability_label_matrix(prob_mat, cutoff_fn = prob_quantile_cutoff)
```

## Arguments

- prob_mat:

  Numeric probability matrix (cells x cell types).

- cutoff_fn:

  Function taking a numeric vector and returning a single numeric cutoff
  (default \[prob_quantile_cutoff()\]).

## Value

Logical matrix with the same dimensions and dimnames as \`prob_mat\`.

## Examples

``` r
prob_mat <- rbind(
  cell1 = c(Bcell = 0.9, Tcell = 0.1),
  cell2 = c(Bcell = 0.2, Tcell = 0.8),
  cell3 = c(Bcell = 0.4, Tcell = 0.3)
)
probability_label_matrix(
  prob_mat,
  cutoff_fn = function(x) 0.5
)
#>       Bcell Tcell
#> cell1  TRUE FALSE
#> cell2 FALSE  TRUE
#> cell3 FALSE FALSE
```
