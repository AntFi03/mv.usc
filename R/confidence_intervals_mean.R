#' Confidence intervals and simultaneous confidence intervals for the mean vector
#'
#' @author Antón Figueroa Martínez
#'
#' @description
#' Function that computes both confidence intervals for summaries of the mean vector given by a*mu, and simultaneous confidence intervals given either by Scheffé or Bonferroni methods.
#'
#' @param dat matrix with the original data the test is based on. Default: `NULL`.
#' @param x.bar sample mean vector. Just in case "dat = `NULL`". Default: `NULL`.
#' @param Sigma_c sample mean corrected matrix. Just in case "dat = `NULL`". Default: `NULL`.
#' @param n number of observacions. Just in case "dat = `NULL`". Default: `NULL`.
#' @param a vector givin the linear combination of mu elements to which compute the confidence interval. Just in case "simultaneous = `NULL`". Default: `NULL`.
#' @param alpha significance level for the interval. Default: `0.05`.
#' @param simultaneous character that selects the method to be applied in the simultaneous confidence intervals computation. Options are: `scheffe` and `bonferroni`. Default: `NULL`.
#' @param verbose bool that enables console summary representation. Default: `TRUE`.
#'
#' @return Returns a matrix with the extremes of the interval(s).
#'
#' @importFrom stats qt
#' @importFrom stats qf
#'
#' @export
confidence_intervals_mean <- function(
  dat = NULL,
  x.bar = NULL,
  Sigma_c = NULL,
  n = NULL,
  a = NULL,
  alpha = 0.05,
  simultaneous = NULL,
  verbose = TRUE
) {
  # -------------------- Preliminaries ------------------- #
  # Check if enough arguments:
  if (missing(dat)) {
    stopifnot(
      "vector x.bar must be specified if dat is missing" = !missing(x.bar),
      "matrix Sigma_c must be specified if dat is missing" = !missing(Sigma_c),
      "integer n must be specified if dat is missing" = !missing(n)
    )
    d <- length(x.bar)
  } else {
    n <- dim(dat)[1]
    d <- dim(dat)[2]
    x.bar <- colMeans(dat)
    Sigma_c <- cov(dat)
  }

  # ------------------------ Cases ----------------------- #
  # Predefine auxiliar function, since it's used in various
  # cases.
  normal_interval <- function(a, x.bar, Sigma_c, n, alpha) {
    # T-value computation:
    tvalue <- qt(alpha / 2, n - 1)

    # Extremes computation:
    ext1 <- t(a) %*% x.bar - tvalue * sqrt((t(a) %*% Sigma_c %*% a) / n)
    ext2 <- t(a) %*% x.bar + tvalue * sqrt((t(a) %*% Sigma_c %*% a) / n)

    # Output:
    output <- list(
      "inf" = min(ext1, ext2),
      "sup" = max(ext1, ext2)
    )

    # End:
    return(output)
  }

  if (missing(simultaneous)) {
    stopifnot(
      "vector a must be specified if simultaneous is missing" = !missing(a),
      "vector a doesn't have the suitable length" = (length(a) == d)
    )

    # --------------- Single Interval Case --------------- #
    # Compute interval:
    output <- normal_interval(a, x.bar, Sigma_c, n, alpha)

    # Transform result:
    result <- list_to_matrix_results(output, rowname = "Mean interval")

    # Message:
    if (verbose) {
      message(
        "\nMultivariate Confidence Interval for the Mean\n"
      )
      message("Analysis results:")
      stats::printCoefmat(
        result,
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
    return(result)
  } else {
    if (!missing(a)) {
      warning(
        "vector a argument will not be used when simultaneous is not null"
      )
    }

    if (simultaneous == "scheffe" || simultaneous == "scheffé") {
      scheffe_interval <- function(
        a,
        x.bar,
        Sigma_c,
        n,
        d,
        alpha
      ) {
        # F-value computation:
        fvalue <- qf(alpha, d, n - d)

        # Extremes computation:
        ext1 <- t(a) %*%
          x.bar -
          sqrt(d * (n - 1) / (n - d) * fvalue) *
            sqrt((t(a) %*% Sigma_c %*% a) / n)
        ext2 <- t(a) %*%
          x.bar +
          sqrt(d * (n - 1) / (n - d) * fvalue) *
            sqrt((t(a) %*% Sigma_c %*% a) / n)

        # Output:
        output <- list(
          "inf" = min(ext1, ext2),
          "sup" = max(ext1, ext2)
        )

        # End:
        return(output)
      }

      # Simultaneous intervals computation:
      output = list()
      for (i in 1:d) {
        # Choose base vectors:
        a <- diag(d)[, i]

        # Compute interval:
        output[[i]] <- scheffe_interval(a, x.bar, Sigma_c, n, d, alpha)
      }

      # Transform results:
      results <- matrix(NA, nrow = 1, ncol = 2)
      for (i in 1:d) {
        results <- rbind(
          results,
          list_to_matrix_results(
            output[[i]],
            rowname = sprintf("X_%d", i)
          )
        )
      }
      results <- results[-1, , drop = FALSE]

      # Message:
      if (verbose) {
        message(
          "\nMultivariate Confidence Simultaneous Scheffé Intervals for the Mean\n"
        )
        message("Analysis results:")
        stats::printCoefmat(
          results,
          digits = 5
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
    } else if (simultaneous == "bonferroni") {
      # Simultaneous intervals computation:
      output = list()
      for (i in 1:d) {
        # Choose base vectors:
        a <- diag(d)[, i]

        # Compute interval:
        output[[i]] <- normal_interval(a, x.bar, Sigma_c, n, alpha / d)
      }

      # Transform results:
      results <- matrix(NA, nrow = 1, ncol = 2)
      for (i in 1:d) {
        results <- rbind(
          results,
          list_to_matrix_results(
            output[[i]],
            rowname = sprintf("X_%d", i)
          )
        )
      }
      results <- results[-1, , drop = FALSE]

      # Message:
      if (verbose) {
        message(
          "\nMultivariate Confidence Simultaneous Bonferroni Intervals for the Mean\n"
        )
        message("Analysis results:")
        stats::printCoefmat(
          results,
          digits = 5
        )
        message("\nMore info:")
        message(sprintf(
          "Number of observations: %d,  Dimension: %d.",
          n,
          d
        ))
      }

      # End:
      return(results)
    } else {
      stop(
        "character simultaneous must take values in {'scheffe', 'bonferroni'}"
      )
    }
  }
}
