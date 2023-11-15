#' Filter a data.table with angles
#'
#' @description
#' Filters out insignificant angles from an angle data.table.
#'
#' @details
#' This is a **low-level** function that utilises functionality
#' of **data.table** package to filter insgificant entries from
#' an angle data.table.
#'
#' @import data.table
#' @param x_df_ang data.table. A long data.table with three columns:
#'   x - name of a gene in row, y - name of a gene in column,
#'   angle - an angle between x and y.
#' @param l_angles list. An output of **estimate.critical.angles**
#'   function.
#' @return list. List with two elements: data.tables with sharp
#' and blunt significant angles.
filter.angles <- function(x_df_ang, #nolint
                          l_angles) {
  ## filter angular data.table by critical angles
  df_ang_sharp <- as.data.table(
    x_df_ang[angle <= l_angles$critical_angles[1]]
    )
  df_ang_sharp <- df_ang_sharp[, edge := paste0(x, "_", y)]
  df_ang_sharp <- df_ang_sharp[, .(edge, angle)]
  #---
  df_ang_blunt <- as.data.table(
    x_df_ang[angle >= l_angles$critical_angles[2]]
  )
  df_ang_blunt <- df_ang_blunt[, edge := paste0(x, "_", y)]
  df_ang_blunt <- df_ang_blunt[, .(edge, angle)]
  invisible(gc())
  ##
  l_critical <- list()
  l_critical[["sharp"]] <- df_ang_sharp
  l_critical[["blunt"]] <- df_ang_blunt
  ##
  return(l_critical)
}