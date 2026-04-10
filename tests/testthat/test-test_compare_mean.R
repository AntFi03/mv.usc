test_that("lecture notes example 1", {
  actual <- test_compare_mean(
    x1.bar = c(136, 102.58, 51.96),
    x2.bar = c(113.38, 88.29, 40.71),
    Sigma1_c = matrix(
      c(432.58, 259.87, 161.67, 259.87, 164.57, 98.99, 161.67, 98.99, 63.87),
      nrow = 3,
      ncol = 3
    ),
    Sigma2_c = matrix(
      c(132.99, 75.85, 35.82, 75.85, 47.96, 20.75, 35.82, 20.75, 10.79),
      nrow = 3,
      ncol = 3
    ),
    n1 = 24,
    n2 = 24
  )

  expect_equal(as.vector(actual$hotelling_statistic), 69.75702, tolerance = 1e-5)
  expect_equal(as.vector(actual$statistic), 22.24137, tolerance = 1e-5)
  expect_equal(as.vector(actual$pvalue), 6.45224e-09, tolerance = 1e-5)
})