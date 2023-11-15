#' Assmble one-hot encoded matrix
#'
#' @description
#' Constructs a matrix where each row is a pair of genes with a significant
#' angle and a column is a sample. Each cell in a matrix can be either 1 or
#' 0, where 1 means that the significant angle between the pair of genes in
#' a row is present in this sample and 0 means not.
#'
#' @import data.table
#' @importFrom magrittr %>%
#' @importFrom purrr map2 reduce
#' @importFrom stringr str_detect
#' @param cons_nodes list. An output from the **assemble.cons.nodes** function.
#' @return data.frame. Rows are gene pairs with significant angles and columns
#' are samples.
#' @export
assemble.onehotmat <- function(cons_nodes) { #nolint
  to_keep <- stringr::str_detect(names(cons_nodes), "conserved", negate = TRUE)
  purrr::map2(
    cons_nodes[to_keep],
    names(cons_nodes[to_keep]),
    function(x, y) {
      tmp <- quote(y)
      tmp <- x[, eval(tmp) := 1]
      return(tmp)
    }
  ) %>%
    purrr::reduce(
      .,
      function(x, y) {
        merge(x, y, all = TRUE, by = "edge")
      }
    ) -> ohm
  ohm[is.na(ohm)] <- 0
  return(ohm)
}