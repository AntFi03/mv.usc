#' Mean vector comparison test between two populations
#'
#' @author Antón Figueroa Martínez
#'
#' @description
#' Function that performs the comparison test between two population mean vectors always when covariance matrices equal and unknown.
#'
#' @param dat1 matrix with the original data the test is based on for the first population. Default: `NULL`.
#' @param dat2 matrix with the original data the test is based on for the second population. Default: `NULL`.
#' @param x1.bar sample mean vector for the first population. Just in case "dat1 = `NULL`". Default: `NULL`.
#' @param x2.bar sample mean vector for the second population. Just in case "dat2 = `NULL`". Default: `NULL`.
#' @param Sigma1_c sample mean corrected matrix for the second population. Just in case "dat1 = `NULL`". Default: `NULL`.
#' @param Sigma2_c sample mean corrected matrix for the second population. Just in case "dat2 = `NULL`". Default: `NULL`.
#' @param n1 number of observacions from the sample of the first population. Just in case "dat1 = `NULL`". Default: `NULL`.
#' @param n2 number of observacions from the sample of the second population. Just in case "dat2 = `NULL`". Default: `NULL`.
#'
#' @return Returns a list with the values of the statistic(s) and p-value.
#'
#' @importFrom stats pf
#'
#' @export
test_compare_mean <- function(
  dat1 = NULL,
  dat2 = NULL,
  x1.bar = NULL,
  x2.bar = NULL,
  Sigma1_c = NULL,
  Sigma2_c = NULL,
  n1 = NULL,
  n2 = NULL
) {
  # -------------------- Preliminaries ------------------- #
  # Check that information is correctly given by pairs:
  stopifnot(
    "data matrices dat1 and dat2 must be both specified or bot missing" = identical(missing(
      dat1
    ),
      missing(dat2)),
    "vectors x1.bar and x1.bar must be both specified or bot missing" = identical(missing(
      x1.bar
    ),
      missing(x2.bar)),
    "matrices Sigma1_c and Sigma2_c must be both specified or bot missing" = identical(missing(
      Sigma1_c
    ),
      missing(Sigma2_c)),
    "parameters n1 and n2 must be both specified or bot missing" = identical(missing(
      n1
    ),
      missing(n2))
  )

  # ------------------------ Cases ----------------------- #
  if (!is.null(dat1)) {
    # ------------------ Case With Data ------------------ #
    # Check if dimensions are comfortable:
    stopifnot(
      "dimensions must be comfortable in dat1 and dat2" = (dim(dat1)[2] ==
        dim(dat2)[2]),
      "dimension should be >1" = (dim(dat1)[2] > 1)
    )

    # Set parameters:
    d <- ncol(dat1)
    n1 <- nrow(dat1)
    n2 <- nrow(dat2)
    x1.bar <- colMeans(dat1)
    x2.bar <- colMeans(dat2)
    Sigma1_c <- cov(dat1)
    Sigma2_c <- cov(dat2)
  } else {
    # ----------------- Case Without Data ---------------- #
    # Check if all arguments are present:
    stopifnot(
      "When non of datX are present, arguments needed: x1.bar, x2.bar, Sigma1_c, Sigma2_c, n1 and n2" = (!missing(
        x1.bar
      ) &&
        !missing(x2.bar) &&
        !missing(Sigma1_c) &&
        !missing(Sigma2_c) &&
        !missing(n1) &&
        !missing(n2))
    )

    # Check if dimensions are comfortable:
    stopifnot(
      "dimensions must be comfortable in x1.bar and x2.bar" = (length(x1.bar) ==
        length(x2.bar)),
      "dimensions must be comfortable in Sigma1_c and Sigma2_c" = ((dim(
        Sigma1_c
      )[
        1
      ] ==
        dim(Sigma1_c)[2]) &&
        (dim(Sigma1_c)[2] == dim(Sigma2_c)[2]) &&
        (dim(Sigma2_c)[1] == dim(Sigma2_c)[2])),
      "dimension should be >1" = (length(x1.bar) > 1)
    )

    # Define d:
    d <- length(x1.bar)
  }

  # ----------------- Common Computations ---------------- #
  # Compute general sample covariance matrix:
  Sigma_c <- ((n1 - 1) * Sigma1_c + (n2 - 1) * Sigma2_c) / (n1 + n2 - 2)

  # Compute Hotelling statistic:
  statistic_hotelling <- ((n1 * n2) / (n1 + n2)) *
    t(x1.bar - x2.bar) %*% solve(Sigma_c) %*% (x1.bar - x2.bar)

  # Compute F statistic:
  statistic_F <- ((n1 + n2 - d - 1) / ((n1 + n2 - 2) * d)) *
    statistic_hotelling

  # Calculate p-value:
  df1 <- d
  df2 <- n1 + n2 - d - 1
  pvalue <- 1 - pf(statistic_F, df1, df2)

  # Output:
  output <- list(
    "hotelling_statistic" = statistic_hotelling,
    "statistic" = statistic_F,
    "pvalue" = pvalue
  )

  # Message:
  message(
    "\nMean Vector Comparison Test Between Two Populations with Unknown Covariance Matrix\n"
  )
  message("Test definition:")
  message("  H0: \u03bc_1 = \u03bc_1")
  message("  Ha: \u03bc_2 != \u03bc_2\n")
  message("Analysis results:")
  stats::printCoefmat(
    list_to_matrix_results(output),
    digits = 5,
    signif.stars = TRUE,
    has.Pvalue = TRUE
  )
  message("\nMore info:")
  message(sprintf(
    "Number of observations: %d and %d,  Dimension: %d,\nF-statistic degrees of freedom: %d and %d.",
    n1,
    n2,
    d,
    df1,
    df2
  ))

  # End:
  return(output)
}
