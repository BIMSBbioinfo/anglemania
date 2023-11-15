#' Read in Seacells from a .RDS file
#'
#' @description
#' import.sl returns a Seurat list with normalised and
#' scaled expression matrices.
#'
#' @details
#' This function doesn't assess if the read object is a
#' Seurat object or not, and therefore relies on the user
#' to provide the path to correct data. Although the function will
#' read and pre-process any Seurat list, metacells are expected as
#' an input.
#'
#' @param path_dat string. Must be a path to an existing
#'   .RDS file with a SeuratList of samples.
#' @param min_mcells integer. Minimal number of Metacells
#'   ppresent in a sample. If less, a sample is discarded.
#'   Defaulit is 15.
#' @return Seurat list.
#' @export
import.sl <- function(path_dat = NULL, min_mcells = 15) { #nolint
  if (is.null(path_dat)) {
    stop("Empty input...\nProvide path to .RDS")
  }
  sl <- readr::read_rds(path_dat)
  # remove unrepresentative datasets
  ds_keep <- purrr::map_dbl(names(sl),
    ~ length(sl[[.x]]$seacell_id)
  ) > min_mcells
  sl <- sl[ds_keep]
  # normalise & scale expression matrix
  sl <- lapply(sl, FUN = NormalizeData)
  sl <- lapply(sl, FUN = ScaleData)
  return(sl)
}