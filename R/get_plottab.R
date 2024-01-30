#' Construct a plottiugn table from a seurat object
#'
#' @description
#' Exctracts a data.frame with embeddings for a specified reduction method
#' and appends specified variables from the meta.data slot.
#'
#' @import Seurat
#' @importFrom magrittr %>%
#' @importFrom dplyr select all_of
#' @param seurat_object seurat object.
#' @param var character vector. A vector of varaibles in meta.data slot
#'   to append.
#' @param reduction character. Name of the reduction method to pull
#'   embeddings from.
#' @return data.frame.
#' @export get_plottab
get_plottab <- function(seurat_object, #nolint
                        var,
                        reduction) {
  cbind(
    seurat_object@meta.data %>% dplyr::select(orig.ident, dplyr::all_of(var)),
    seurat_object@reductions[[reduction]]@cell.embeddings
  )
}