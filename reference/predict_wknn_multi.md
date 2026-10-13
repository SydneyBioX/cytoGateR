# Multi-Metric Weighted kNN

Multi-Metric Weighted kNN

## Usage

``` r
predict_wknn_multi(
  train_data,
  test_data,
  train_labels,
  k = 5,
  method = "pearson",
  return_matrix = FALSE,
  chunk_size = 10000L
)
```

## Arguments

- train_data:

  Data frame or matrix of training cells (rows = markers, cols = cells).

- test_data:

  Data frame or matrix of test cells (rows = markers, cols = cells).

- train_labels:

  Factor of labels for the training data.

- k:

  Number of neighbors.

- method:

  One of "pearson", "spearman", "cosine", or "euclidean"

- return_matrix:

  Logical indicating whether to return the full probability matrix for
  all classes. If `FALSE`, only the predicted labels and associated
  probabilities are returned. Default is `FALSE`.

- chunk_size:

  Optional integer number of test cells per block. If `NULL`, all test
  cells are processed in a single block (legacy behavior).

## Value

A list containing predicted 'labels' and 'probs'.
