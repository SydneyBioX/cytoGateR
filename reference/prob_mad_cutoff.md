# Median + MAD cutoff helper

Computes a cutoff as median(x) + \`mad_mult\` \* MAD(x).

## Usage

``` r
prob_mad_cutoff(x, mad_mult = 3)
```

## Arguments

- x:

  Numeric vector of probabilities.

- mad_mult:

  Multiplier on MAD (default 3).

## Value

Numeric scalar cutoff.

## Examples

``` r
prob_mad_cutoff(c(0.1, 0.1, 0.2, 0.9), mad_mult = 2)
#> [1] 0.25
```
