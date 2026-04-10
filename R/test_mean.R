#' Mean multivariate test
#'
#' @author Antón Figueroa Martínez
#'
#'
#' @description
#' Function that contrast a simple hypothesis over the means for multivariate normal populations when the covariance matrix is know and when it is not. It can also handle generalized mean tests based on A*mu=b hypothesis.
#'
#'
#' @param dat data.frame or matrix to which apply the test. Default: `NULL`.
#' @param mu0 mean vector under null hypothesis. Default: `NULL`, but mandatory if not generalized test.
#' @param cov.matrix population covariance matrix (if known). If not specified, substituted by the sample covariance matrix (with the corresponding adjustment in the test design). Default: `NULL`.'
#' @param n number of observations. Just in case "dat = `NULL`". Default: `NULL`.'
#' @param x.bar sample mean vector. Just in case "dat = `NULL`". Default: `NULL`.
#' @param A matrix giving the restriction coefficients for a generalized test. Default: `NULL`.
#' @param b vector giving the right-sides of the restrictions in a generalized test. Default: `NULL`.
#' @param verbose bool that enables console summary representation. Default: `TRUE`.
#'
#'
#' @return Returns the value of the statistic and the p-value of the test.
#'
#' @importFrom stats pchisq
#' @importFrom stats pf
#' @importFrom stats cov
#' @importFrom stats printCoefmat
#'
#' @export
#'
#' @examples
#' # Data:
#' Sigma <- matrix(c(49, 35, 35, 52), nrow = 2)
#' x.bar <- c(102, 98)
#' mu0 <- c(100, 100)
#'
#' # Function:
#' test_mean(n = 20, cov.matrix = Sigma, x.bar = x.bar, mu0 = mu0)
#'
test_mean <- function(
  dat = NULL,
  mu0 = NULL,
  cov.matrix = NULL,
  n = NULL,
  x.bar = NULL,
  A = NULL,
  b = NULL,
  verbose = TRUE
) {
  # -------------------- Preliminaries ------------------- #
  # Auto-set parameters when data is given:
  if (!is.null(dat)) {
    n <- nrow(dat)
    x.bar <- colMeans(dat)
  } else {
    stopifnot(
      "matrix cov.matrix must be specified when no data is available" = !is.null(
        cov.matrix
      ),
      "parameter n must be specified when no data is available" = !is.null(n),
      "vector x.bar must be specified when no data is available" = !is.null(
        x.bar
      )
    )
  }

  stopifnot(
    "vector b must be specified to perform a generalized test (matrix A alone is not enough" = !(!missing(
      A
    ) &&
      missing(b)),
    "matrix A must be specified to perform a generalized test (vector b alone is not enough" = !(!missing(
      b
    ) &&
      missing(A))
  )

  # ------------------------ Cases ----------------------- #
  # Is the contrast a common one or a generalized one
  # (matrix A and vector b of restrictions must be given):
  if (is.null(A) & is.null(b)) {
    # Stop if mu0 missing (nothing to compare with):
    stopifnot(
      "vector m0 must be specified when the test is not a generalized one (without it, there is nothing to compare to)" = !missing(
        mu0
      )
    )

    # Set dimension:
    d <- length(mu0)

    # Is cov.matrix specified?:
    if (!is.null(cov.matrix)) {
      # --------------- Case With Known Sigma -------------- #
      Sigma <- cov.matrix

      # Statistic computation:
      statistic <- n * (x.bar - mu0) %*% solve(Sigma) %*% (x.bar - mu0)

      # P-value computation:
      df <- d
      pvalue <- 1 - pchisq(statistic, df)

      # Output:
      output <- list("statistic" = statistic, "pvalue" = pvalue)

      # Message:
      if (verbose) {
        message(
          "\nMean Simple Hypothesis Test with Known Covariance Matrix\n"
        )
        message("Test definition:")
        message("  H0: \u03bc = \u03bc_0")
        message("  Ha: \u03bc != \u03bc_0\n")
        message("Analysis results:")
        stats::printCoefmat(
          list_to_matrix_results(output),
          digits = 5,
          signif.stars = TRUE,
          has.Pvalue = TRUE
        )
        message("\nMore info:")
        message(sprintf(
          "Number of observations: %d,  Dimension: %d,\nChi-statistic degrees of freedom: %d.",
          n,
          d,
          df
        ))
      }

      # End:
      return(output)
    } else {
      # -------------- Case With Unknown Sigma ------------- #
      stopifnot(
        "matrix dat is needed since the original data is used to compute sample covariance matrix" = !missing(
          dat
        )
      )
      Sigma_c <- cov(dat)
      Sigma_e <- (n - 1) / n * Sigma_c

      # Hotelling statistic computation:
      statistic_hotelling <- n *
        (x.bar - mu0) %*% solve(Sigma_c) %*% (x.bar - mu0)

      # F statistic computation:
      statistic_F <- (n - d) /
        d *
        (x.bar - mu0) %*% solve(Sigma_e) %*% (x.bar - mu0)

      # P-value computation:
      df1 <- d
      df2 <- n - d
      pvalue <- 1 - pf(statistic_F, df1 = df1, df2 = df2)

      # Output:
      output <- list(
        "hotelling_statistic" = statistic_hotelling,
        "statistic" = statistic_F,
        "pvalue" = pvalue
      )

      # Message:
      if (verbose) {
        message(
          "\nMean Simple Hypothesis Test with Unknown Covariance Matrix\n"
        )
        message("Test definition:")
        message("  H0: \u03bc = \u03bc_0")
        message("  Ha: \u03bc != \u03bc_0\n")
        message("Analysis results:")
        stats::printCoefmat(
          list_to_matrix_results(output),
          digits = 5,
          signif.stars = TRUE,
          has.Pvalue = TRUE
        )
        message("\nMore info:")
        message(sprintf(
          "Number of observations: %d,  Dimension: %d,\nF-statistic degrees of freedom: %d and %d.",
          n,
          d,
          df1,
          df2
        ))
      }

      # End:
      return(output)
    }
    # Is the test a generalized one?
  } else if (!is.null(A) & !is.null(b)) {
    # ----------------- Generalized Test ----------------- #
    # Warn about unused arguments:
    if (!missing(mu0)) {
      warning(
        "vector mu0 is not being used, since the function is performing a generalized test"
      )
    }

    # Set dimension:
    if (!is.null(Sigma)) {
      d <- dim(Sigma)[2]
    } else {
      stopifnot(
        "matrix dat is needed since the original data is used to compute sample covariance matrix" = !missing(
          dat
        )
      )
      d <- dim(dat)[2]
    }

    # Check matrix A and vector b dimensions:
    stopifnot(
      "uncomfortable dimensions for matrix A and vector b" = !((dim(as.matrix(
        A
      ))[1] ==
        length(b)) &&
        (dim(as.matrix(A))[2] == d))
    )
    q <- dim(A)[1]

    # Is cov.matrix specified?:
    if (!is.null(cov.matrix)) {
      # ------- (Generalized) Case With Known Sigma ------ #
      Sigma <- cov.matrix

      # Statistic computation:
      statistic <- n *
        t(A %*% x.bar - b) %*% solve(A %*% Sigma %*% t(A)) %*% (A %*% x.bar - b)

      # P-value computation:
      df <- q
      pvalue <- 1 - pchisq(statistic, df)

      # Output:
      output <- list("statistic" = statistic, "pvalue" = pvalue)

      # Message:
      if (verbose) {
        message(
          "\nMean Generalized Test with Known Covariance Matrix\n"
        )
        message("Test definition:")
        message("  H0: A\u03bc = b")
        message("  Ha: A\u03bc != b\n")
        message("Analysis results:")
        stats::printCoefmat(
          list_to_matrix_results(output),
          digits = 5,
          signif.stars = TRUE,
          has.Pvalue = TRUE
        )
        message("\nMore info:")
        message(sprintf(
          "Number of observations: %d,  Dimension: %d,\nChi-statistic degrees of freedom: %d.",
          n,
          d,
          df
        ))
      }

      # End:
      return(output)
    } else {
      # ------ (Generalized) Case With Unknown Sigma ----- #
      Sigma_c <- cov(dat)
      Sigma_e <- (n - 1) / n * Sigma_c

      # Hotelling statistic computation:
      statistic_hotelling <- (n - 1) *
        t(A %*% x.bar - b) %*%
          solve(A %*% Sigma_e %*% t(A)) %*%
          (A %*% x.bar - b)

      # F statistic computation:
      statistic_F <- (n - q) /
        q *
        t(A %*% x.bar - b) %*%
          solve(A %*% Sigma_e %*% t(A)) %*%
          (A %*% x.bar - b)

      # P-value computation:
      df1 <- q
      df2 <- n - q
      pvalue <- 1 - pf(statistic_F, df1 = df1, df2 = df2)

      # Output:
      output <- list(
        "hotelling_statistic" = statistic_hotelling,
        "statistic" = statistic_F,
        "pvalue" = pvalue
      )

      # Message:
      if (verbose) {
        message(
          "\nMean Generalized Test with Unknown Covariance Matrix\n"
        )
        message("Test definition:")
        message("  H0: A\u03bc = b")
        message("  Ha: A\u03bc != b\n")
        message("Analysis results:")
        stats::printCoefmat(
          list_to_matrix_results(output),
          digits = 5,
          signif.stars = TRUE,
          has.Pvalue = TRUE
        )
        message("\nMore info:")
        message(sprintf(
          "Number of observations: %d,  Dimension: %d,\nF-statistic degrees of freedom: %d and %d.",
          n,
          d,
          df1,
          df2
        ))
      }

      # End:
      return(output)
    }
  }
}
