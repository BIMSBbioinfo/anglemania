#' Integrate Samples in a Seurat Object Using Selected Features from \code{anglem_object} Object
#'
#' @description
#' `integrate_by_features` integrates samples or batches within a Seurat object using canonical correlation analysis (CCA) based on a set of selected features (genes). The function utilizes an `anglem` object to extract integration genes and handles the integration process, including optional downstream processing steps such as scaling, PCA, and UMAP visualization.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item **Batch Key Addition**: Adds a unique batch key to the Seurat object's metadata to distinguish different batches or samples. Batch key is set to the \code{anglem_object}'s \code{batch_key}.
#'   \item **Splitting**: Splits the Seurat object into a list of Seurat objects based on the batch key.
#'   \item **Integration**: Calls \code{\link{integrate_seurat_list}} to integrate the list of Seurat objects using the features extracted from the \code{anglem_object}.
#' }
#'
#' The integration is performed using Seurat's CCA-based methods, and parameters are adjusted based on the smallest dataset to ensure compatibility with small sample sizes (e.g., metacells or SEACells). If \code{process = TRUE}, the function will also scale the data, run PCA, and compute UMAP embeddings.
#'
#' @param seurat_object A \code{\link[Seurat]{Seurat}} object containing all samples or batches to be integrated.
#' @param anglem_object An \code{\link{anglem}} object previously generated using \code{\link{create_anglem}} and \code{\link{big_anglemanise}}. It is important that the \code{dataset_key} and \code{batch_key} are correctly set in the \code{anglem_object}.
#' @param int_order An optional data frame specifying the integration order of samples within the Seurat list. See the \code{sample.tree} argument in \code{\link[Seurat]{IntegrateData}} for more details. If not provided, Seurat will construct the integration order using hierarchical clustering. Default is \code{NULL}.
#' @param process Logical value indicating whether to further process the data after integration (i.e., scale it, run PCA, and compute UMAP embeddings). Default is \code{TRUE}.
#' @param verbose Logical value indicating whether to display progress messages during integration. Default is \code{FALSE}.
#'
#' @return A \code{\link[Seurat]{Seurat}} object containing the integrated data. The default assay is set to \code{"integrated"}.
#'
#' @importFrom Seurat SplitObject NormalizeData FindIntegrationAnchors IntegrateData ScaleData RunPCA RunUMAP DefaultAssay
#' @importFrom pbapply pblapply
#'
#' @seealso
#' \code{\link{create_anglem}},
#' \code{\link{big_anglemanise}},
#' \code{\link{extract_integration_genes}},
#' \code{\link{integrate_seurat_list}},
#' \code{\link[Seurat]{IntegrateData}},
#' \code{\link[Seurat]{FindIntegrationAnchors}}
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(Seurat)
#' library(pbapply)
#'
#' # Assume seurat_object is a Seurat object containing all samples/batches
#' # and anglem_object is an anglem object with integration genes identified
#'
#' # Integrate the Seurat object using the anglem object
#' integrated_seurat <- integrate_by_features(
#'   seurat_object = seurat_object,
#'   anglem_object = anglem_object,
#'   process = TRUE,
#'   verbose = TRUE
#' )
#'
#' # Access the integrated assay
#' DefaultAssay(integrated_seurat) # Should be "integrated"
#'
#' # Visualize UMAP embedding
#' DimPlot(integrated_seurat, reduction = "umap", group.by = "batch")
#' }
#'
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


#' Integrate a List of Seurat Objects Using Selected Features
#'
#' @description
#' `integrate_seurat_list` integrates a list of Seurat objects (e.g., representing different samples or batches) using canonical correlation analysis (CCA) based on a set of selected features (genes). The function handles normalization, finding integration anchors, integrating data, and optional downstream processing steps such as scaling, PCA, and UMAP visualization.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item **Normalization**: Each Seurat object in the list is log-normalized using \code{\link[Seurat]{NormalizeData}}.
#'   \item **Parameter Adjustment**: Integration parameters are adjusted based on the smallest dataset to accommodate cases with a small number of cells (e.g., metacells).
#'   \item **Finding Integration Anchors**: Uses \code{\link[Seurat]{FindIntegrationAnchors}} to find anchors between datasets based on the provided features.
#'   \item **Integration**: Integrates the datasets using \code{\link[Seurat]{IntegrateData}}.
#'   \item **Optional Processing**: If \code{process = TRUE}, the function scales the data, runs PCA, and computes UMAP embeddings.
#' }
#'
#' The integration is performed using Seurat's CCA-based methods, and the function is designed to handle datasets with varying sizes efficiently.
#'
#' @param seurat_list A list of \code{\link[Seurat]{Seurat}} objects to be integrated.
#' @param features A character vector of gene names (features) used for integration.
#' @param int_order An optional data frame specifying the integration order of samples within the Seurat list. See the \code{sample.tree} argument in \code{\link[Seurat]{IntegrateData}} for more details. If not provided, Seurat will construct the integration order using hierarchical clustering. Default is \code{NULL}.
#' @param process Logical value indicating whether to further process the data after integration (i.e., scale it, run PCA, and compute UMAP embeddings). Default is \code{TRUE}.
#' @param verbose Logical value indicating whether to display progress messages during integration. Default is \code{FALSE}.
#'
#' @return A \code{\link[Seurat]{Seurat}} object containing the integrated data. The default assay is set to \code{"integrated"}.
#'
#' @importFrom Seurat NormalizeData FindIntegrationAnchors IntegrateData ScaleData RunPCA RunUMAP DefaultAssay
#' @importFrom pbapply pblapply
#'
#' @seealso
#' \code{\link{integrate_by_features}},
#' \code{\link[Seurat]{IntegrateData}},
#' \code{\link[Seurat]{FindIntegrationAnchors}},
#' \code{\link[Seurat]{NormalizeData}},
#' \code{\link[Seurat]{ScaleData}},
#' \code{\link[Seurat]{RunPCA}},
#' \code{\link[Seurat]{RunUMAP}}
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(Seurat)
#' library(pbapply)
#'
#' # Assume seurat_list is a list of Seurat objects representing different samples
#' # and features is a vector of gene names used for integration
#'
#' # Integrate the Seurat list
#' integrated_seurat <- integrate_seurat_list(
#'   seurat_list = seurat_list,
#'   features = features,
#'   process = TRUE,
#'   verbose = TRUE
#' )
#'
#' # Access the integrated assay
#' DefaultAssay(integrated_seurat) # Should be "integrated"
#'
#' # Visualize UMAP embedding
#' DimPlot(integrated_seurat, reduction = "umap", group.by = "batch")
#' }
#'
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
