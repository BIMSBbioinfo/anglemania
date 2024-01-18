#' Write angles as tab-delimted tables
#'
#' @description
#' Writes to the disk two tables with angles in a tab-delimited format.
#'
#' @details
#' This is a **low-level** function that uses **data.table's** package
#' fwrite function to write tables with croitical angles to the disk.
#' Names of the written tables are md5-sums of the same tables.
#'
#' @importFrom data.table fwrite
#' @importFrom digest digest
#' @param l_critical list. List with two elements: data.tables with sharp
#' and blunt significant angles.
#' @param path_to_write_angles string. Path to a directory to write tables.
#' @return list. List with md5-sums correspondign to the two written tables.
#' @export write_angles
write_angles <- function(l_critical, #nolint
                         path_to_write_angles) {
  ##
  l_df_ang_md5 <- list()
  for (i in names(l_critical)) {
    md5 <- digest::digest(l_critical[[i]], "md5")
    data.table::fwrite(
      x = l_critical[[i]],
      file = file.path(path_to_write_angles, md5),
      sep = "\t"
    )
    l_df_ang_md5[[i]] <- md5
  }
  return(l_df_ang_md5)
}