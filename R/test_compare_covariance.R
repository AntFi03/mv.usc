#' Covariance matrix comparison test between two populations
#'
#' @author Antón Figueroa Martínez
#'
#' @description
#' Function that performs the comparison test between k normal population covariance matrices.
#'
#' @param dat_list list of matrices with the original data for each population's sample. Default: `NULL`.
#' @param Sigma_c_list list of sample mean corrected matrices for the each population. Just in case "dat_list = `NULL`". Default: `NULL`.
#' @param n_list list of numbers of observacions from the sample of each population. Just in case "dat_list = `NULL`". Default: `NULL`.
#' @param verbose bool that enables console summary representation. Default: `TRUE`.
#'
#' @return Returns a list with the values of the statistic and p-value.
#'
#' @importFrom stats pchisq
#'
#' @export
test_compare_covariance <- function(
  dat_list = NULL,
  Sigma_c_list = NULL,
  n_list = NULL,
  verbose = TRUE
) {
  # -------------------- Preliminaries ------------------- #

  # ------------------------ Cases ----------------------- #
  if (!is.null(dat_list)) {
    # ------------------ Case With Data ------------------ #
    # Get list length:
    k <- length(dat_list)
    d <- ncol(dat_list[[1]])

    # Check if dimensions are comfortable:
    stopifnot(
      "dimensions must be comfortable in dat_list matrices" = (length(unique(sapply(
        dat_list,
        ncol
      ))) ==
        1),
      "dimension should be >1" = (dim(dat_list[[1]])[2] > 1)
    )

    # Set parameters:
    n_list <- vapply(dat_list, nrow, integer(1))
    Sigma_c_list <- lapply(dat_list, cov)
  } else {
    # ----------------- Case Without Data ---------------- #
    # Check if all arguments are present:
    stopifnot(
      "When dat_list is missing, arguments needed: Sigma_c_list and n_list" = (!missing(
        Sigma_c_list
      ) &&
        !missing(n_list))
    )

    # Check that information is correctly given by lists with
    # the same length:
    stopifnot(
      "length of all lists must be the same" = identical(
        length(Sigma_c_list),
        length(n_list)
      )
    )

    # Get list length:
    k <- length(n_list)
    d <- ncol(Sigma_c_list[[1]])

    # Check if dimensions are comfortable:
    stopifnot(
      "dimensions must be comfortable in Sigma_c_list" = all(vapply(
        Sigma_c_list,
        function(x) {
          identical(dim(x), dim(Sigma_c_list[[1]])) && dim(x)[1] == dim(x)[2]
        },
        logical(1)
      )),
      "dimension should be >1" = (dim(Sigma_c_list[[1]])[2] > 1)
    )
  }

  # ----------------- Common Computations ---------------- #
  # Compute general sample covariance matrix:
  Q_matrix <- Reduce(`+`, Map(`*`, Sigma_c_list, n_list))
  n <- Reduce(`+`, n_list)
  Qn_matrix <- Q_matrix / n

  # Compute the statistic:
  statistic <- sum(vapply(
    seq_along(Sigma_c_list),
    function(i) {
      n_list[[i]] *
        as.numeric(
          determinant(Qn_matrix, logarithm = TRUE)$modulus -
            determinant(Sigma_c_list[[i]], logarithm = TRUE)$modulus
        )
    },
    numeric(1)
  ))

  # Calculate p-value:
  df <- 0.5 * d * (d + 1) * (k - 1)
  pvalue <- 1 - pchisq(statistic, df)

  # Output:
  output <- list(
    "statistic" = statistic,
    "pvalue" = pvalue
  )

  # Message:
  if (verbose) {
    message(
      "\nCovariance Matrix Comparison Test Between k Normal Populations\n"
    )
    message("Test definition:")
    message("  H0: \u03A3_1 = ... = \u03A3_k")
    message("  Ha: \u03A3_i != \u03A3_j\n")
    message("Analysis results:")
    stats::printCoefmat(
      list_to_matrix_results(output),
      digits = 5,
      signif.stars = TRUE,
      has.Pvalue = TRUE
    )
    message("\nMore info:")
    message(sprintf(
      "Mean number of observations: %d,  Dimension: %d,\nChi-statistic degrees of freedom: %d.",
      n / k,
      d,
      df
    ))
  }

  # End:
  return(output)
}
