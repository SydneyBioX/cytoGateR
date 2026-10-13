# Compute the cutoff where posteriors are equal

Compute the cutoff where posteriors are equal

## Usage

``` r
gmm_equal_posterior_cutoff(mu1, mu2, s1, s2, p1, p2)
```

## Arguments

- mu1:

  Numeric mean of component 1.

- mu2:

  Numeric mean of component 2.

- s1:

  Numeric SD of component 1.

- s2:

  Numeric SD of component 2.

- p1:

  Numeric mixing weight of component 1.

- p2:

  Numeric mixing weight of component 2.

## Value

Numeric cutoff. Defaults to midpoint if root finding fails.
