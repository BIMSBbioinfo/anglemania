#' Assemble conserved matrix
#'
#' @description
#' Constructs a matrix where each row is a pair of genes with a significant
#' angle and a column is a sample. Each cell in a matrix can be either 1 or
#' 0, where 1 means that the significant angle between the pair of genes in
#' a row is present in this sample and 0 means not.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfr
#' @param ohm data.frame. An output from the **assemble_onehotmat** function.
#' @return data.frame. Rows are gene pairs with significant angles and columns
#' are samples.
#' @export assemble_conmat
assemble_conmat <- function(ohm) { #nolint
  idcs <- colnames(ohm[, 2:dim(ohm)[2]])
  ord <- order(colSums(ohm[, ..idcs]), decreasing = TRUE)
  idcs <- idcs[ord]
  purrr::map_dfr(
    idcs,
    function(idx) {
      tmp <- ohm[, ..idcs]
      tmp <- tmp[tmp[[idx]] > 0, ]
      out <- colSums(tmp)
      return(out)
    }
  ) %>%
    as.matrix(.) -> m_ass
  rownames(m_ass) <- colnames(m_ass)
  return(m_ass)
}