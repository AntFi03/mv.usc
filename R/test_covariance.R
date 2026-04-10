#' Covariance multivariate tests
#'
#' @author Antón Figueroa Martínez
#'
#' @description
#' Function that performs some contrasts related to the covariance matrix of a normal multivariate sample.
#'
#' @param dat data.frame or matrix to which apply the test. Default: `NULL`.
#' @param Sigma_0 covariance matrix under null hypothesis. Mandatory for "simple" and "proportional" tests. Default:`NULL`
#' @param Sigma_e Sample covariance matrix. Default: `NULL`.
#' @param n number of observations. Just in case "dat = `NULL`". Default: `NULL`.
#' @param tests Vector to determine the tests to perform. Could take integer and character values according to: {1: "simple", 2: "proportional", 3: "independence"}.
#' @param verbose bool that enables console summary representation. Default: `TRUE`.
#'
#' @return Returns a list with the statistics and p-values of the corresponding tests.
#'
#' @importFrom stats pchisq
#' @importFrom stats pf
#' @importFrom stats printCoefmat
#' @importFrom stats cov2cor
#'
#' @export
test_covariance <- function(
  dat = NULL,
  Sigma_0 = NULL,
  Sigma_e = NULL,
  n = NULL,
  tests = c("simple"),
  verbose = TRUE
) {
  # -------------------- Preliminaries ------------------- #
  # Auto-set parameters when data is given:
  if (!is.null(dat)) {
    n <- nrow(dat)
    Sigma_c <- cov(dat)
    Sigma_e <- (n - 1) / n * Sigma_c

    # Set dimension:
    d <- ncol(as.matrix(dat))
    stopifnot(
      "dimension d must be >1, since makes no sense to test covariance matrix when there is not one" = (d >
        1)
    )
  } else {
    stopifnot(
      "parameter n must be specified when no data is available" = !is.null(n),
      "matrix Sigma_e must be specified when no data is available" = !is.null(
        Sigma_e
      )
    )

    # Set dimension:
    d <- ncol(as.matrix(Sigma_e))
    stopifnot(
      "dimension d must be >1, since makes no sense to test covariance matrix when there is not one" = (d >
        1)
    )

    Sigma_c <- n / (n - 1) * Sigma_e
  }

  # Case determination:
  possible_tests <- c(
    "simple",
    "proportional",
    # "subvector_independence",
    "independence"
  )
  possible_tests_long <- c(
    "Covariance Matrix Simple Hypothesis Equality Test",
    "Covariance Matrix Simple Hypothesis Proportionality Test",
    # "Covariance Matrix Subvector Independence Test",
    "Covariance Matrix Independence Test"
  )

  if ("all" %in% tests) {
    chosen_tests <- rep(TRUE, 3)
  } else {
    chosen_tests <- sapply(seq_along(possible_tests), function(i) {
      i %in% tests | possible_tests[i] %in% tests
    })
  }

  output = list()

  # --------------------- Five Cases --------------------- #
  if (chosen_tests[1]) {
    case = 1
    stopifnot(
      "matrix argument Sigma_0 must be specified for simple test" = !missing(
        Sigma_0
      )
    )
    # -------- Simple Covariance Matrix Comparison ------- #
    lambda <- eigen(solve(Sigma_0) %*% Sigma_e)$values
    a <- mean(lambda)
    g <- prod(lambda)^(1 / d)

    # Statistic computation:
    statistic <- n * d * (a - log(g) - 1)

    # P-value computation:
    df <- d * (d + 1) / 2
    pvalue <- 1 - pchisq(statistic, df = df)

    # Output:
    output[[case]] <- list(
      "statistic" = statistic,
      "df" = df,
      "pvalue" = pvalue
    )
    names(output)[case] <- possible_tests[case]
  }
  if (chosen_tests[2]) {
    case = 2
    stopifnot(
      "matrix argument Sigma_0 must be specified for simple test" = !missing(
        Sigma_0
      )
    )
    # ------ Covariance Matrix Proportionality Test ------ #
    lambda <- eigen(solve(Sigma_0) %*% Sigma_e)$values
    a0 <- mean(lambda)
    g0 <- prod(lambda)^(1 / d)

    # Statistic computation:
    statistic <- n * d * log(a0 / g0)

    # P-value computation:
    df <- (d - 1) * (d + 2) / 2
    pvalue <- 1 - pchisq(statistic, df = df)

    # Output:
    output[[case]] <- list(
      "statistic" = statistic,
      "df" = df,
      "pvalue" = pvalue
    )
    names(output)[case] <- possible_tests[case]
  }
  if (chosen_tests[3]) {
    case = 3
    # -------- Covariance Matrix Independency Test ------- #
    # Statistic computation:
    statistic <- -n * log(det(cov2cor(Sigma_c)))

    # P-value computation:
    df <- d * (d - 1) / 2
    pvalue <- 1 - pchisq(statistic, df = df)

    # Output:
    output[[case]] <- list(
      "statistic" = statistic,
      "df" = df,
      "pvalue" = pvalue
    )
    names(output)[case] <- possible_tests[case]
  }

  if (sum(chosen_tests) == 1) {
    # Get rid of list (because just one output):
    names(output) <- NULL
    output <- unlist(output, recursive = FALSE)

    # Message:
    if (verbose) {
      message(
        sprintf("\n%s\n", possible_tests_long[chosen_tests])
      )
      message("Test definition:")
      switch(
        possible_tests[chosen_tests],
        "simple" = message(
          "  H0: \u03A3 = \u03A3_0\n  Ha: \u03A3 != \u03A3_0\n"
        ),
        "proportional" = message(
          "  H0: \u03A3 = k*\u03A3_0\n  Ha: \u03A3 != k*\u03A3_0\n"
        ),
        # "subvector_independence" = message(
        #   "  H0: \u03A3_{12} = 0\n  Ha: \u03A3_{12} != 0\n"
        # ),
        "independence" = message(
          "  H0: \u03A3 diagonal\n  Ha: \u03A3 free\n"
        )
      )
      message("Analysis results:")
      stats::printCoefmat(
        list_to_matrix_results(
          output,
          rowname = sprintf("Test_%s", possible_tests[chosen_tests])
        ),
        digits = 5,
        signif.stars = TRUE,
        has.Pvalue = TRUE
      )
      message("\nMore info:")
      message(sprintf(
        "Number of observations: %d,  Dimension: %d.",
        n,
        d
      ))
    }

    # End:
    return(output)
  } else {
    # Message:
    if (verbose) {
      message(
        sprintf("\nSet of Covariance Matrix Multivariate Tests:")
      )
      results <- matrix(NA, nrow = 1, ncol = 3)
      for (i in seq_along(chosen_tests)) {
        if (chosen_tests[i]) {
          message(
            sprintf("  - %s", possible_tests_long[i])
          )
          results <- rbind(
            results,
            list_to_matrix_results(
              output[[i]],
              rowname = sprintf("Test %s", possible_tests[i])
            )
          )
        }
      }
      results <- results[-1, , drop = FALSE]
      # results <- results[1:end]
      # message("\nTest definition:")
      # switch(
      #   possible_tests[chosen_tests],
      #   "simple" = message(
      #     "  H0: \u03A3 = \u03A3_0\n  Ha: \u03A3 != \u03A3_0\n"
      #   ),
      #   "proportional" = message(
      #     "  H0: \u03A3 = k*\u03A3_0\n  Ha: \u03A3 != k*\u03A3_0\n"
      #   ),
      #   "subvector independence" = message(
      #     "  H0: \u03A3_{12} = 0\n  Ha: \u03A3_{12} != 0\n"
      #   ),
      #   "independence" = message(
      #     "  H0: \u03A3 diagonal\n  Ha: \u03A3 free\n"
      #   )
      # )
      message("\nAnalysis results:")
      stats::printCoefmat(
        results,
        digits = 5,
        signif.stars = TRUE,
        has.Pvalue = TRUE
      )
      message("\nMore info:")
      message(sprintf(
        "Number of observations: %d,  Dimension: %d.",
        n,
        d
      ))
    }

    # End:
    return(output)
  }
}
