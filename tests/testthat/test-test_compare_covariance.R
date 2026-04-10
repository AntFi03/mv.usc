test_that("lecture notes example 1", {
  actual <- test_compare_covariance(
  Sigma_c_list = list(
    matrix(
      c(432.58, 259.87, 161.67, 259.87, 164.57, 98.99, 161.67, 98.99, 63.87),
      nrow = 3,
      ncol = 3
    ),
    matrix(
      c(132.99, 75.85, 35.82, 75.85, 47.96, 20.75, 35.82, 20.75, 10.79),
      nrow = 3,
      ncol = 3
    )
  ),
  n_list = list(24, 24)
)

  expect_equal(as.vector(actual$statistic), 26.98611, tolerance = 1e-5)
  expect_equal(as.vector(actual$pvalue), 0.0001456775, tolerance = 1e-5)
})