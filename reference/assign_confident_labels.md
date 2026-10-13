# Assign labels with uncertainty handling

Use a probability matrix and confidence flags to assign lineage labels;
cells that are not confident receive \`unknown_label\`.

## Usage

``` r
assign_confident_labels(prob_mat, confidence_flags, unknown_label = "Unknown")
```

## Arguments

- prob_mat:

  Numeric matrix of probabilities (cells x lineages).

- confidence_flags:

  Logical vector length \`nrow(prob_mat)\`.

- unknown_label:

  Character label for uncertain cells (default "Unknown").

## Value

Character vector of labels length \`nrow(prob_mat)\`.

## Examples

``` r
prob_mat <- rbind(
  cell1 = c(Bcell = 0.9, Tcell = 0.1),
  cell2 = c(Bcell = 0.5, Tcell = 0.5)
)
assign_confident_labels(prob_mat, c(TRUE, FALSE))
#> [1] "Bcell"   "Unknown"
```
