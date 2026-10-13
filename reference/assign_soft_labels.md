# Assign hard labels from probability matrix

Assign hard labels from probability matrix

## Usage

``` r
assign_soft_labels(prob_mat, unknown_thresh = 0.4)
```

## Arguments

- prob_mat:

  Matrix \[cells x types\].

- unknown_thresh:

  If max prob \< this, label as "Unknown" (default 0.4).

## Value

Character vector of labels length nrow(prob_mat).

## Examples

``` r
prob_mat <- rbind(
  cell1 = c(Bcell = 0.8, Tcell = 0.2),
  cell2 = c(Bcell = 0.3, Tcell = 0.35)
)
assign_soft_labels(prob_mat, unknown_thresh = 0.4)
#> [1] "Bcell"   "Unknown"
```
