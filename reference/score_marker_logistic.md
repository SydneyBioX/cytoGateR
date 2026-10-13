# Logistic marker score

Logistic marker score

## Usage

``` r
score_marker_logistic(x, cutoff, scale)
```

## Arguments

- x:

  numeric

- cutoff:

  numeric

- scale:

  numeric

## Value

numeric in \[0,1\]

## Examples

``` r
score_marker_logistic(c(0, 1, 2), cutoff = 1, scale = 0.5)
#> [1] 0.1192029 0.5000000 0.8807971
```
