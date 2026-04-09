#' Auxiliary reshape results function
#' 
#' @author Antón Figueroa Martínez (antfigtel@gmail.com)
#' 
#' @description
#' Function takes the results of a test (list type) to matrix type for results printing.
#' 
#' @param results_list List with test results.
#' @param rowname Name for the test to be identified with. Default = "Test".
#' 
#' @return Returns a matrix with the ready-to-print results.
#'  
list_to_matrix_results <- function(results_list, rowname = "Test") {
  # List to row matrix transformation:
  results_row <- do.call(cbind, results_list)
  colnames(results_row) <- tools::toTitleCase(names(results_list))
  colnames(results_row)[colnames(results_row) == "Pvalue"] <- "Pr(>|t|)"
  rownames(results_row) <- rowname

  # End:
  return(results_row)
}