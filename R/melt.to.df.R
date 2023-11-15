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
#' @import data.table
#' @param x_mat_ang matrix. Contains angles between the
#' vectors of gene expression between genes.
#' @return data.table. A long data.table with three columns:
#' x - name of a gene in row, y - name of a gene in column,
#' angle - an angle between x and y.
melt.to.df <- function(x_mat_ang) { #nolint
  ##
  x_mat_ang[upper.tri(x_mat_ang, diag = TRUE)] <- NA
  x_df_ang <- data.table::as.data.table(x_mat_ang, keep.rownames = "x")
  x_df_ang <- data.table::melt(x_df_ang,
                               id.vars = "x",
                               variable.name = "y",
                               value.name    = "angle")
  x_df_ang <- na.omit(x_df_ang)
  x_df_ang <- x_df_ang[, y := as.character(y)]
  class(x_df_ang) <- "data.table"
  return(x_df_ang)
}