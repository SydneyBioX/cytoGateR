test_that("marker scoring functions return expected values", {
  expect_equal(
    score_marker_logistic(c(0, 1, 2), cutoff = 1, scale = 0.5),
    stats::plogis(c(-2, 0, 2))
  )
  expect_equal(score_marker_rank(c(3, 1, 2)), c(1, 0, 0.5))
})

test_that("path scores and tree probability follow the selected branch", {
  tree <- list(
    type = "node",
    marker = "CD3e",
    cutoff = 1,
    scale = 0.5,
    left = list(type = "leaf"),
    right = list(type = "leaf")
  )
  expr_mat <- rbind(CD3e = c(cell1 = 0.5, cell2 = 1.5))

  expect_equal(
    collect_path_scores(tree, expr_mat, cell_i = 1),
    stats::plogis(-1)
  )
  expect_equal(
    tree_prob(tree, expr_mat, cell_i = 2, lambda = 0),
    stats::plogis(1)
  )
})

test_that("negative marker expression decreases a score", {
  expr_mat <- rbind(CD3e = c(low = 0, high = 2))
  marker_stats <- list(CD3e = list(cutoff = 1, scale = 0.5))

  low_penalty <- neg_penalty(
    expr_mat, "CD3e", 1,
    marker_stats = marker_stats,
    neg_strength = 0.8
  )
  high_penalty <- neg_penalty(
    expr_mat, "CD3e", 2,
    marker_stats = marker_stats,
    neg_strength = 0.8
  )

  expect_gt(low_penalty, high_penalty)
  expect_true(all(c(low_penalty, high_penalty) >= 0))
  expect_true(all(c(low_penalty, high_penalty) <= 1))
})
