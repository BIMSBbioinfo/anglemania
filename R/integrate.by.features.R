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
#' @return seurat object. Integrated seurat object. Seurat::DefaultAssay
#' is set to "integrated".
#' @export
integrate.by.features <- function(seurat_list,
                                  features_to_integrate,
                                  int_order = NULL) {
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = seurat_list,
    anchor.features = features_to_integrate,
    verbose = FALSE,
    dims = 1:10,
    k.filter  = 10,
    k.anchor  = 10,
    k.score   = 10,
    reduction = "cca"
  )
  features_intersect <- Reduce(
    function(x, y) intersect(x, y), lapply(seurat_list, rownames)
  )
  seurat_combined <- Seurat::IntegrateData(
    anchorset = anchors,
    k.weight = 10,
    features.to.integrate = features_intersect,
    sample.tree = int_order,
    verbose = TRUE
  )
  Seurat::DefaultAssay(seurat_combined) <- "integrated"
  seurat_combined <- Seurat::ScaleData(seurat_combined)
  seurat_combined <- Seurat::RunPCA(seurat_combined, npcs = 30)
  seurat_combined <- Seurat::RunUMAP(
    seurat_combined,
    reduction = "pca",
    dims = 1:30
  )
  return(seurat_combined)
}