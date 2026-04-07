#' Mean multivariate test
#'
#' @author Antón Figueroa Martínez
#' (antfigtel@gmail.com)
#'
#'
#' @description Function that contrast a simple hypothesis over the means for multivariate normal populations when the covariance matrix is know and when it is not.
#'
#'
#' @param dat data.frame or matriz to which apply the test. Default: `NULL`.
#'
#' @param mu0 mean vector under null hypothesis. Mandatory
#'
#' @param cov.matrix population covariance matrix (if known). If not specified, substituted by the sample covariance matrix (with the corresponding adjustment in the test design). Default: `NULL`.
#'
#' @param alpha pre-specified significance level. Default: 0.05.
#'
#' @param n number of observations. Just in case "dat = `NULL`". Default: `NULL`.
#'
#' @param x.bar sample mean vector. Just in case "dat = `NULL`". Default: `NULL`.
#'
#'
#' @return Returns the value of the statistic and the p-value of the test.
#'
#' @importFrom stats pchisq
#' @importFrom stats pf
#' @importFrom stats cov
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
  mu0,
  cov.matrix = NULL,
  alpha = 0.05,
  n = NULL,
  x.bar = NULL
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
      "vector x.bar must be specified when no data is available" = !is.null(x.bar)
    )
  }

  # Set dimension:
  d <- length(mu0)

  # ---------------------- Two Cases --------------------- #
  # Is cov.matrix specified?:
  if (!is.null(cov.matrix)) {
    # --------------- Case With Known Sigma -------------- #
    Sigma <- cov.matrix

    # Statistic computation:
    statistic <- n * (x.bar - mu0) %*% solve(Sigma) %*% (x.bar - mu0)

    # P-value computation:
    pvalue <- 1 - pchisq(statistic, d)

    # Output:
    return(list("statistic" = statistic, "pvalue" = pvalue))
  } else {
    # -------------- Case With Unknown Sigma ------------- #
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
    pvalue <- 1 - pf(statistic_F, df1 = d, df2 = n - d)

    # Output:
    return(list(
      "statistic" = statistic_F,
      "pvalue" = pvalue,
      "hotelling_statistic" = statistic_hotelling
    ))
  }
}
