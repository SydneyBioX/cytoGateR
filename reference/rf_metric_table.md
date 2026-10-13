# Tabular and text summaries of RF performance

Convenience helpers to expose overall metrics, cross-validation
summaries, and per-class statistics in both tabular and text-friendly
formats.

## Usage

``` r
rf_metric_table(fit)
```

## Arguments

- fit:

  List returned by \[train_custom_rf()\].

## Value

\`rf_metric_table()\` returns a list of tables: \`overall\`,
\`cv_overall\`, \`cv_overall_summary\`, \`class\`, and \`class_raw\`.
\`rf_metric_text()\` returns a character vector of summary lines.

## Examples

``` r
fit <- list(
  metrics = list(
    accuracy = 0.75,
    f1_macro = 0.73,
    cv_overall = data.frame(
      accuracy = c(0.7, 0.8),
      f1_macro = c(0.68, 0.78)
    )
  ),
  test_truth = c("Bcell", "Bcell", "Tcell", "Tcell"),
  test_pred = c("Bcell", "Tcell", "Tcell", "Tcell")
)
rf_metric_table(fit)
#> $overall
#>     metric value
#> 1 accuracy  0.75
#> 2 f1_macro  0.73
#> 
#> $cv_overall
#>   accuracy f1_macro
#> 1      0.7     0.68
#> 2      0.8     0.78
#> 
#> $cv_overall_summary
#>   accuracy_mean accuracy_sd f1_macro_mean f1_macro_sd
#> 1          0.75  0.07071068          0.73  0.07071068
#> 
#> $class
#>       class tp fp fn support precision recall        f1
#> Bcell Bcell  1  0  1       2 1.0000000    0.5 0.6666667
#> Tcell Tcell  2  1  0       2 0.6666667    1.0 0.8000000
#> 
#> $class_raw
#>       class tp fp fn support precision recall        f1
#> Bcell Bcell  1  0  1       2 1.0000000    0.5 0.6666667
#> Tcell Tcell  2  1  0       2 0.6666667    1.0 0.8000000
#> 
```
