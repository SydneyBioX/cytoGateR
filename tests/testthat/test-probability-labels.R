test_that("confidence flags and labels are assigned consistently", {
  prob_mat <- rbind(
    cell1 = c(Bcell = 0.90, Tcell = 0.10),
    cell2 = c(Bcell = 0.55, Tcell = 0.45),
    cell3 = c(Bcell = 0.05, Tcell = 0.95)
  )

  flags <- compute_custom_labels(
    prob_mat,
    flag_fn = function(x) max(x) >= 0.8
  )

  expect_identical(flags, c(TRUE, FALSE, TRUE))
  expect_identical(
    assign_confident_labels(prob_mat, flags),
    c("Bcell", "Unknown", "Tcell")
  )
})

test_that("soft labels respect the unknown threshold", {
  prob_mat <- rbind(
    cell1 = c(Bcell = 0.8, Tcell = 0.2),
    cell2 = c(Bcell = 0.3, Tcell = 0.35)
  )

  expect_identical(
    assign_soft_labels(prob_mat, unknown_thresh = 0.4),
    c("Bcell", "Unknown")
  )
})

test_that("probability cutoffs produce a logical label matrix", {
  prob_mat <- rbind(
    cell1 = c(Bcell = 0.9, Tcell = 0.1),
    cell2 = c(Bcell = 0.2, Tcell = 0.8),
    cell3 = c(Bcell = 0.4, Tcell = 0.3)
  )

  label_mat <- probability_label_matrix(
    prob_mat,
    cutoff_fn = function(x) 0.5
  )

  expect_type(label_mat, "logical")
  expect_identical(dimnames(label_mat), dimnames(prob_mat))
  expect_identical(rowSums(label_mat), c(cell1 = 1, cell2 = 1, cell3 = 0))
  expect_equal(prob_quantile_cutoff(1:4, prob = 0.5), 2.5)
  expect_true(is.numeric(prob_mad_cutoff(c(0.1, 0.2, 0.3))))
})
