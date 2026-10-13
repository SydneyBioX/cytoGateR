# Flag confident cells from a probability matrix

Compute per-cell confidence using a user-supplied decision function or a
robust Tukey-style rule: a cell is considered confident if its maximum
lineage probability exceeds \`max(base_thresh, quantile + iqr_mult \*
IQR)\` for that cell.

## Usage

``` r
compute_custom_labels(
  prob_mat,
  base_thresh = 0.4,
  quantile_cut = 0.75,
  iqr_mult = 1.5,
  flag_fn = NULL
)
```

## Arguments

- prob_mat:

  Numeric matrix of probabilities (cells x lineages).

- base_thresh:

  Minimum allowable threshold (default 0.4) used only for the default
  rule.

- quantile_cut:

  Quantile used to set the adaptive threshold (default 0.75) for the
  default rule.

- iqr_mult:

  Multiplier applied to the IQR (default 1.5) for the default rule.

- flag_fn:

  Optional function taking a numeric vector of probabilities and
  returning a single logical value. Defaults to the Tukey-style rule
  above.

## Value

Logical vector length \`nrow(prob_mat)\` indicating confident cells.

## Examples

``` r
prob_mat <- rbind(
  cell1 = c(Bcell = 0.9, Tcell = 0.1),
  cell2 = c(Bcell = 0.5, Tcell = 0.5)
)
compute_custom_labels(
  prob_mat,
  flag_fn = function(x) max(x) >= 0.8
)
#> [1]  TRUE FALSE
```
