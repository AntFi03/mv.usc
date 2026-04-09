test_that("lecture notes example 1", {
  actual <- test_covariance(
    dat = read.table(
      "http://eamo.usc.es/pub/pateiro/files/AM_learnr/data//notas.txt",
      header = TRUE
    ),
    Sigma_0 = matrix(c(49, 35, 35, 52), nrow = 2),
    # Sigma_e = S,
    # n = n,
    test = c(1)
  )

  expect_equal(as.vector(actual$statistic), 10.47213, tolerance = 1e-5)
  expect_equal(as.vector(actual$df), 3, tolerance = 1e-5)
  expect_equal(as.vector(actual$pvalue), 0.01495, tolerance = 1e-4)
})