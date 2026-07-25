test_that("class metrics calculate expected confusion counts", {
  truth <- c("Bcell", "Bcell", "Tcell", "Tcell")
  pred <- c("Bcell", "Tcell", "Tcell", "Tcell")
  metrics <- class_metrics_from_fit(
    list(),
    test_truth = truth,
    test_pred = pred
  )

  bcell <- metrics[as.character(metrics$class) == "Bcell", ]
  tcell <- metrics[as.character(metrics$class) == "Tcell", ]

  expect_equal(bcell$tp, 1)
  expect_equal(bcell$fn, 1)
  expect_equal(tcell$tp, 2)
  expect_equal(tcell$fp, 1)
})

test_that("F1 calculation returns one row per reference class", {
  data("cytoGateR_example", package = "cytoGateR")
  result <- calculate_f1(
    cytoGateR_example,
    ref_col = "cell_type_hard",
    pred_col = "knn_label"
  )

  expected_types <- unique(as.character(
    SummarizedExperiment::colData(cytoGateR_example)$cell_type_hard
  ))
  expect_setequal(as.character(result$Category), expected_types)
  expect_true(all(result$F1_Score >= 0 & result$F1_Score <= 1))
})
