# Derive per-class metrics from a fitted model

Prefers cross-validation class metrics when available; otherwise derives
metrics from the provided truth and prediction vectors.

## Usage

``` r
class_metrics_from_fit(fit, test_truth = NULL, test_pred = NULL)
```

## Arguments

- fit:

  List returned by \[train_custom_rf()\].

- test_truth:

  Optional truth vector (used only if CV metrics absent).

- test_pred:

  Optional prediction vector (used only if CV metrics absent).

## Value

Data frame with columns \`class\`, \`tp\`, \`fp\`, \`fn\`, \`support\`,
\`precision\`, \`recall\`, \`f1\`, and optionally \`fold\` when sourced
from CV.

## Examples

``` r
truth <- c("Bcell", "Bcell", "Tcell", "Tcell")
pred <- c("Bcell", "Tcell", "Tcell", "Tcell")
class_metrics_from_fit(list(), truth, pred)
#>       class tp fp fn support precision recall        f1
#> Bcell Bcell  1  0  1       2 1.0000000    0.5 0.6666667
#> Tcell Tcell  2  1  0       2 0.6666667    1.0 0.8000000
```
