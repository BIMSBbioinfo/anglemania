#' Construct a plottiugn table from a seurat object
#'
#' @description
#' Exctracts a data.frame with embeddings for a specified reduction method
#' and appends specified variables from the meta.data slot.
#'
#' @import Seurat
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of
#' @param seu_obj seurat object.
#' @param var character vector. A vector of varaibles in meta.data slot
#'   to append.
#' @param reduction character. Name of the reduction method to pull
#'   embeddings from.
#' @return data.frame.
#' @export
get.plottab <- function(seu_obj, #nolint
                        var,
                        reduction) {
  cbind(
    seu_obj@meta.data %>% dplyr::select(orig.ident, dplyr::all_of(var)),
    seu_obj@reductions[[reduction]]@cell.embeddings
  )
}