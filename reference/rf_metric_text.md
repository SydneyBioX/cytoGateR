# Text summary of random forest performance

Formats key overall and cross-validation metrics from a fitted model
returned by \[train_custom_rf()\] into a character vector of
human-readable lines.

## Usage

``` r
rf_metric_text(fit, digits = 3)
```

## Arguments

- fit:

  List returned by \[train_custom_rf()\].

- digits:

  Number of digits when formatting numeric values (default 3).

## Value

A character vector of summary lines.

## Examples

``` r
fit <- list(
  metrics = list(
    accuracy = 0.75,
    f1_macro = 0.73,
    cv_overall = NULL
  ),
  test_truth = c("Bcell", "Bcell", "Tcell", "Tcell"),
  test_pred = c("Bcell", "Tcell", "Tcell", "Tcell")
)
rf_metric_text(fit)
#> [1] "Overall accuracy: 0.750"   "Overall F1 (macro): 0.730"
```
