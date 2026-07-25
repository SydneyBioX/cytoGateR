test_that("packaged example data contain required analysis outputs", {
  data("cytoGateR_example", package = "cytoGateR")

  expect_s4_class(cytoGateR_example, "SpatialExperiment")
  expect_true("exprs" %in% SummarizedExperiment::assayNames(cytoGateR_example))
  expect_equal(nrow(cytoGateR_example), 40)
  expect_gt(ncol(cytoGateR_example), 0)

  cd_names <- names(SummarizedExperiment::colData(cytoGateR_example))
  expect_true(all(c(
    "sample_id", "image_name", "core_group",
    "knn_label", "knn_confidence"
  ) %in% cd_names))
})

test_that("stored kNN probabilities are complete and normalized", {
  data("cytoGateR_example", package = "cytoGateR")
  cd <- SummarizedExperiment::colData(cytoGateR_example)
  prob_cols <- grep("^KNN_P_", names(cd), value = TRUE)
  prob_mat <- as.matrix(cd[, prob_cols])

  expect_length(prob_cols, 10)
  expect_false(anyNA(prob_mat))
  expect_equal(
    unname(rowSums(prob_mat)),
    rep(1, nrow(prob_mat)),
    tolerance = 1e-8
  )
})

test_that("uncertainty output is aligned to cells", {
  data("cytoGateR_example", package = "cytoGateR")
  spe <- cytoGateR_example[, seq_len(12L)]
  cd <- SummarizedExperiment::colData(spe)
  prob_cols <- grep("^KNN_P_", names(cd), value = TRUE)
  prob_mat <- as.matrix(cd[, prob_cols])
  rownames(prob_mat) <- colnames(spe)
  colnames(prob_mat) <- sub("^KNN_P_", "", prob_cols)

  uncertainty <- calculate_uncertainty(
    prob_mat,
    spe,
    sample_col = "sample_id",
    k_spatial = 20
  )

  expect_identical(uncertainty$cell_id, colnames(spe))
  expect_equal(nrow(uncertainty), ncol(spe))
  expect_true(all(uncertainty$entropy >= 0))
  expect_true(all(uncertainty$entropy <= 1))
})
