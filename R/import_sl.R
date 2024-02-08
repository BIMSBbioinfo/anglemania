#' Read in Seacells from a .RDS file
#'
#' @description
#' import_sl takes a list of Seurat objects, checks if the conditions for 
#' anglemanise are correct.
#' Briefly, it checks 
#' a) if all objects have more than 'mcells'. Discard those samples with less.
#' b) a layer with scaled data. If not, it performs SCTransform (using glampoi)
#' returns a list of corrected Seurat objects
#' @details
#' This function doesn't assess if the read object is a
#' Seurat object or not, and therefore relies on the user
#' to provide the path to correct data. Although the function will
#' read and pre-process any Seurat list, metacells are expected as
#' an input.
#' @importFrom SeuratObject Layers
#' @importFrom Seurat SCTransform
#' @param seurat_list Seurat list.
#' @param min_mcells integer. Minimal number of Metacells
#'   present in a sample. If less, a sample is discarded.
#'   Default is 6. Integration using Seurat CCA seemingly
#'   doesn't work with <6 cells.
#' @return Seurat list.
#' @export import_sl
import_sl <- function(seurat_list, min_mcells = 6) { #nolint
  
  ####### Validate Input ########
  if (length(seurat_list)<2){
    stop("Seurat list must at least contain two Seurat objects")
  }
  if (any(sapply(seurat_list, function(x) attr(class(x), "package") != "SeuratObject"))){
    stop("List contains elements that are not Seurat objects")
  }
  if (!is.integer(min_mcells) || min_mcells <= 5) {
    stop("The argument 'min_mcells' must be an integer larger than 5.")
  }
  
  
  
  ####### PREPARE FOR ANGLEMANISE  ##########
  ## remove unrepresentative datasets

  message(paste0("Removing data sets with less than ", min_mcells, " SEACells"))
  ds_keep <- sapply(seurat_list, function(x) ncol(x) >= min_mcells)
  seurat_list <- seurat_list[ds_keep]
  
  
  ##  check if scale.data is present and SCTransform is necessary
  if (!all(sapply(seurat_list, function(x) "scale.data" %in% SeuratObject::Layers(x)))){
    message("scaled data is not present in all samples. 
  Applying SCTransform to all samples.")
    seurat_list <- lapply(seurat_list,
                          function(x) {
                            Seurat::SCTransform(x,
                                                return.only.var.genes = FALSE,
                                                min_cells = 2,
                                                verbose = FALSE)
                          }
    )
  }
  return(seurat_list)
}
