test_that("lecture notes example 1", {
  actual <- test_mean(
    n = 20,
    cov.matrix = matrix(c(49, 35, 35, 52), nrow = 2),
    x.bar = c(102, 98),
    mu0 = c(100, 100)
  )

  expect_equal(as.vector(actual$statistic), 10.34014, tolerance = 1e-5)
  expect_equal(as.vector(actual$pvalue), 0.005684182, tolerance = 1e-7)
})


test_that("lecture notes example 2", {
  actual <- test_mean(
    dat = read.table(
      paste(
        "http://eamo.usc.es/pub/pateiro/files/AM_learnr/data/",
        "notas.txt",
        sep = ""
      ),
      header = TRUE
    ),
    cov.matrix = matrix(c(49, 35, 35, 52), nrow = 2),
    mu0 = c(100, 100)
  )

  expect_equal(as.vector(actual$statistic), 10.34014, tolerance = 1e-5)
  expect_equal(as.vector(actual$pvalue), 0.005684182, tolerance = 1e-7)
})


test_that("error message check", {
  expect_error(
    test_mean(
      n = 20,
      # cov.matrix = matrix(c(49, 35, 35, 52), nrow = 2),
      x.bar = c(102, 98),
      mu0 = c(100, 100)
    ),
    regexp = "matrix cov.matrix must be specified when no data is available"
  )
})