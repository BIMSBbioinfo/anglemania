#' Extract gene featuries from the processed datasets
#'
#' @description
#' Extract unique gene features per dataset with significant
#' sharp/blunt angles.
#'
#' @importFrom purrr map
#' @param l_processed list. An output from the **anglemanise** function.
#' @return list. List with unique gene features and conserved angles
#' between genes.
#' @export featurise
featurise <- function(l_processed) {
  cons_nodes <- assemble_cons_nodes(
    l_processed,
    c("blunt", "sharp")
  )
  l_ufts <- purrr::map(
    cons_nodes,
    ~ count_unqiue_features(.x$conserved)
  )
  ##
  out <- list()
  out$unique_features <- l_ufts
  out$conserved_nodes <- cons_nodes
  return(l_ufts)
}