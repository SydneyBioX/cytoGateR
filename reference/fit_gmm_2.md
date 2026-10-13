# Fit a 2-component Gaussian mixture model

Fit a 2-component Gaussian mixture model

## Usage

``` r
fit_gmm_2(
  x,
  cutoff_method = c("mean", "equal_posteriors"),
  gmm_model_names = NULL
)
```

## Arguments

- x:

  numeric vector

- cutoff_method:

  Character. One of "mean" (midpoint of means) or "equal_posteriors"
  (threshold where component posteriors are equal).

- gmm_model_names:

  Optional character vector of model names to pass to
  \[mclust::Mclust()\] (e.g., "V" to forbid equal-variance in 1D).

## Value

list with GMM parameters or NULL

## Examples

``` r
set.seed(1)
x <- c(rnorm(50, 0, 0.2), rnorm(50, 2, 0.2))
fit <- fit_gmm_2(x)
fit$cutoff
#> [1] 1.021777
```
