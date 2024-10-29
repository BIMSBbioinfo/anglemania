#' Integrate samples in a seurat list
#'
#' @description
#' Integrates samples in the seurat list using CCA and parameters devised
#' by VF.
#'
#' @import Seurat
#' @param seurat_object seurat object containing all samples/batches.
#' @param anglem_object anglem object. Previously generated using create_anglem() and big_anglemanise(). Important to set the right dataset_key and batch_key!
#' @param features_to_integrate character vector. Vector of gene names
#'   (features) to construct anchors from.
#' @param int_order data.frame. Table specifying the integration order
#'   of samples within the seurat list. See Seurat::IntegrateData's argument
#'   sample.tree help page for more details. If not provided Seurat will
#'   construct the integration order from hclust. Default is NULL.
#' @param process bool. Boolean TRUE/FALSE which indicates whether to further
#'   process the data, i.e., scale it and run PCA & UMAP. Default is TRUE.
#' @return seurat object. Integrated seurat object. Seurat::DefaultAssay
#' is set to "integrated".
#' @export integrate_by_features
integrate_by_features <- function(seurat_object,
                                  anglem_object,
                                  int_order = NULL,
                                  process = TRUE, # should scaling, PCA, and UMAP already be performed on the integrated data?
                                  verbose = FALSE) {

  dataset_key <- anglem_object@dataset_key
  batch_key <- anglem_object@batch_key
  seurat_object <- add_unique_batch_key(seurat_object, dataset_key, batch_key) # temporarily adds unique batch key "batch" to metadata of the seurat object
  seurat_list <- Seurat::SplitObject(seurat_object, split.by = "batch") # split by batch

  seurat_combined <- integrate_seurat_list(
    seurat_list = seurat_list,
    features = extract_integration_genes(anglem_object),
    int_order = int_order,
    process = process,
    verbose = verbose
  )

  return(seurat_combined)
}



#' Integrate samples in a seurat list
#' 
#' @import Seurat
#' @param seurat_list list of seurat objects containing all samples/batches.
#' @param features character vector. Vector of gene names used for integration
#' @param int_order data.frame. Table specifying the integration order
#'   of samples within the seurat list. See Seurat::IntegrateData's argument
#'   sample.tree help page for more details. If not provided Seurat will
#'   construct the integration order from hclust. Default is NULL.
#' @param process bool. Boolean TRUE/FALSE which indicates whether to further
#'   process the data, i.e., scale it and run PCA & UMAP. Default is TRUE.
#' @return seurat object. Integrated seurat object. Seurat::DefaultAssay
#' is set to "integrated".
#' @export integrate_seurat_list
integrate_seurat_list <- function(seurat_list,
                                  features,
                                  int_order = NULL,
                                  process = TRUE, # should scaling, PCA, and UMAP already be performed on the integrated data?
                                  verbose = FALSE) {

  message("Log normalizing data...")
  seurat_list <- pblapply(seurat_list, Seurat::NormalizeData, verbose = verbose)
  smallest_DS <- min(sapply(seurat_list, function(x) ncol(x))) # When working with SEACells/metacells, the number of cells can be really low. Adjust the following algorithm parameters accordingly.
  n <- ifelse(smallest_DS <= 10, smallest_DS - 1, 10)
  message("Finding integration anchors...")
  anchors <- Seurat::FindIntegrationAnchors( # COMMENT: don't worry about normalizing and scaling, the function does log normalization and scaling by default! https://satijalab.org/seurat/reference/findintegrationanchors
    object.list = seurat_list,
    anchor.features = features,
    verbose = verbose,
    dims = 1:n,
    k.filter = n,
    k.anchor = n,
    k.score = n,
    reduction = "cca"
  )

  message("Integrating samples...")
  seurat_combined <- Seurat::IntegrateData( # automatically performs log normalization BUT NOT SCALING
    anchorset = anchors,
    dims = 1:n,
    k.weight = n,
    sample.tree = int_order,
    verbose = verbose
  )

  if (process) {
    get_n_neighbors <- function(num_cells) {
      if (num_cells < 30) {
        return(5)
      } else if (num_cells >= 30 && num_cells <= 100) {
        return(10)
      } else {
        return(30)
      }
    }
    #### RUN PCA AND UMAP ####
    n_seacells <- ncol(seurat_combined)
    n_pcs <- ifelse(n_seacells <= 30, n_seacells - 1, 30)
    n_neighbors <- get_n_neighbors(n_seacells)
    message("Running PCA with ", n_pcs, " PCs")
    
    Seurat::DefaultAssay(seurat_combined) <- "integrated"
    seurat_combined <- Seurat::ScaleData(seurat_combined, verbose = verbose)
    seurat_combined <- Seurat::RunPCA(seurat_combined, npcs = n_pcs, verbose = verbose)

    message("Running UMAP with ", n_pcs, " PCs and ", n_neighbors, " neighbors")
    seurat_combined <- Seurat::RunUMAP(
      seurat_combined,
      reduction = "pca",
      dims = 1:n_pcs,
      n.neighbors = n_neighbors,
      verbose = verbose
    )
  }

  return(seurat_combined)
}
