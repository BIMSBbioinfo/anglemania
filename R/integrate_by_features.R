#' Integrate samples in a seurat list
#'
#' @description
#' Integrates samples in the seurat list using CCA and parameters devised
#' by VF.
#'
#' @import Seurat
#' @param seurat_list seurat list.
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
integrate_by_features <- function(seurat_list,
                                  features_to_integrate,
                                  int_order = NULL,
                                  process = TRUE) {
  smallest_DS <- min(sapply(seurat_list, function(x) ncol(x)))
  # seu_intgr <- integrate_by_features(sl, purrr::reduce(l_ifts, union))    # possible to provide custom integration order
  n <- ifelse(smallest_DS <= 10, smallest_DS - 1, 10)
  print(n)
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = seurat_list,
    anchor.features = features_to_integrate,
    verbose = TRUE,
    dims = 1:n,
    k.filter = n,
    k.anchor = n,
    k.score = n,
    reduction = "cca"
  )

  features_intersect <- Reduce(
    function(x, y) intersect(x, y), lapply(seurat_list, rownames)
  )

  seurat_combined <- Seurat::IntegrateData(
    anchorset = anchors,
    dims = 1:n,
    k.weight = n,
    features.to.integrate = features_intersect,
    sample.tree = int_order,
    verbose = TRUE
  )


  Seurat::DefaultAssay(seurat_combined) <- "integrated"
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
    #### SCALE DATA ####
    seurat_combined <- Seurat::SCTransform(seurat_combined,
      return.only.var.genes = FALSE,
      min_cells = 5, # default
      verbose = FALSE
    )
    #### RUN PCA AND UMAP ####
    n_seacells <- ncol(seurat_combined)
    n_pcs <- ifelse(n_seacells <= 30, n_seacells - 1, 30)
    n_neighbors <- get_n_neighbors(n_seacells)

    message("Running PCA with ", n_pcs, " PCs")
    seurat_combined <- Seurat::RunPCA(seurat_combined, npcs = n_pcs)
    message("Running UMAP with ", n_pcs, " PCs and ", n_neighbors, " neighbors")
    seurat_combined <- Seurat::RunUMAP(
      seurat_combined,
      reduction = "pca",
      dims = 1:n_pcs,
      n.neighbors = n_neighbors
    )
  }


  return(seurat_combined)
}
