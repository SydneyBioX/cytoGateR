# Marker separability

Marker separability

## Usage

``` r
marker_separability(x)
```

## Arguments

- x:

  numeric

## Value

numeric

## Examples

``` r
set.seed(1)
x <- c(rnorm(50, 0, 0.2), rnorm(50, 2, 0.2))
marker_separability(x)
#> [1] 7.92586
```
