# Negative marker penalty

Compute a multiplicative penalty factor in the range \[0, 1\] for a
given cell based on the expression of predefined negative markers.

## Usage

``` r
neg_penalty(
  expr_mat,
  neg_markers,
  cell_i,
  marker_stats = NULL,
  neg_strength = 0.8
)
```

## Arguments

- expr_mat:

  A numeric matrix-like object containing expression values, with
  markers in rows and cells in columns.

- neg_markers:

  A character vector of marker names to penalize.

- cell_i:

  An integer index specifying the column (cell) to evaluate.

- marker_stats:

  A named list containing per-marker statistics. Each element must
  include at least:

  cutoff

  :   Numeric value representing the logistic midpoint.

  scale

  :   Numeric value controlling the logistic slope.

- neg_strength:

  A numeric scalar in \[0, 1\] controlling the strength of the penalty.
  Larger values increase the impact of negative-marker expression.

## Value

A numeric scalar in the range \[0, 1\] representing the multiplicative
penalty factor for the specified cell.

## Details

For each negative marker, a logistic scoring function is used to
estimate the probability (\`p_bad\`) that the marker is expressed above
its expected background level. The final penalty is computed as the
product of (1 - p_bad \* neg_strength) across all negative markers.

Higher negative-marker expression results in a stronger penalty (i.e., a
smaller multiplicative factor).

If no valid negative markers are found in `marker_stats`, the function
returns 1 (no penalty).

## Examples

``` r
expr_mat <- rbind(CD3 = c(0.2, 2))
marker_stats <- list(CD3 = list(cutoff = 1, scale = 0.5))
neg_penalty(
  expr_mat,
  neg_markers = "CD3",
  cell_i = 1,
  marker_stats = marker_stats
)
#> [1] 0.8656147
```
