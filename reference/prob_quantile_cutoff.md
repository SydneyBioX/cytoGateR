# Default quantile-based cutoff

Computes a per-vector cutoff using the specified quantile.

## Usage

``` r
prob_quantile_cutoff(x, prob = 0.98)
```

## Arguments

- x:

  Numeric vector of probabilities.

- prob:

  Quantile to use (default 0.98).

## Value

Numeric scalar cutoff.

## Examples

``` r
prob_quantile_cutoff(c(0.1, 0.2, 0.8, 0.9), prob = 0.75)
#> [1] 0.825
```
