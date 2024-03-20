#' Melt distance matrix
#'
#' @description
#' Melts an angle matrix into a long data.frame.
#'
#' @details
#' This is a **low-level** function that utilises functionality
#' of **data.table** pacckage to melt a square angle matrix into
#' a long data.frame.
#'
#' @importFrom data.table as.data.table melt
#' @param x_mat_ang matrix. Contains angles between the
#' vectors of gene expression between genes.
#' @return data.table. A long data.table with three columns:
#' x - name of a gene in row, y - name of a gene in column,
#' angle - an angle between x and y.
#' @export melt_to_df
melt_to_df <- function(x_mat_ang) { #nolint
  ##
  x_mat_ang[upper.tri(x_mat_ang, diag = TRUE)] <- NA
  
  idx <- which(!is.na(x_mat_ang), arr.ind = TRUE)
  
  # Create a data table from indices and values
  x_df_ang <- data.table::data.table(
    x = rownames(x_mat_ang)[idx[, 1]],
    y = colnames(x_mat_ang)[idx[, 2]],
    angle = x_mat_ang[idx]
  )
  
  return(x_df_ang)
}