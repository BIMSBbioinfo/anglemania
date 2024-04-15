#' Record critical angles between genes in gene expression
#' matrices from a Seurat list.
#'
#' @description
#' *anglemanise* is a high-level function that returns a list
#' of lists containing critical angles between genes across
#' input samples, angle statistics per sample, paths to the
#' the angle records per sample.
#'
#' @details
#' This is the main function of the package. It calculates
#' angles between all genes across all samples provided in
#' the Seurat list. It approximates the angle distribution
#' for each sample and extracts values of critical angles. Data needs to be scaled.
#'
#' @importFrom Matrix Matrix
#' @importFrom purrr map
#' @importFrom SeuratObject LayerData
#' @param seurat_list seurat list of seurat objects with scaled data.
#' @param extrema double. Fraction of the angles
#'   to be cut from both sides of an approximated angle
#'   distribution.
#' @param n_threads integer. Number of cores to use for
#'   the computation of the cosine distances.
#' @param path_to_write_angles string. A path, under which
#'   a temporary directory will be created to record values
#'   of critical angles for each sample.
#' @return list. First two elements are sparse matrices
#'   recording the number of sharp and blunt angles across
#'   the datasets. Third and fourth elements are lists with
#'   angles statistics and paths to the values of critical
#'   angles.
#' @seealso https://arxiv.org/abs/1306.0256
#' @export anglemanise
anglemanise <- function(seurat_list, #nolint
                        extrema = 0.001,
                        n_threads = 16,
                        path_to_write_angles = ".") {
  
  ############## Validate inputs ###########################
  
  if (!is.list(seurat_list) || any(sapply(seurat_list, function(x) attr(class(x), "package") != "SeuratObject"))) {
    stop("seurat_list needs to be a list of Seurat objects")
  }
  if (!is.numeric(n_threads) || n_threads < 1) {
    stop("n_threads has to be a positive integer")
  }
  if (!is.character(path_to_write_angles)) {
    stop("path_to_write_angles has to be a string")
  }
  if (length(extrema) > 2 || length(extrema) < 1) {
    stop("extrema should be numeric of length either 1 or 2")
  }
  ##  check if scale.data is present and SCTransform is necessary
  if (!all(sapply(seurat_list, function(x) "scale.data" %in% SeuratObject::Layers(x)))){
    stop("scaled data is not present in all samples. \nPlease run import_sl on the Seurat list")
  }
  
  if (length(extrema) == 1) {
    extrema <- c(extrema, 1 - extrema)
  }

  if (!dir.exists(path_to_write_angles)) {
    dir.create(path_to_write_angles, recursive = TRUE)
  }
  
  ## intersect genes ==> matrix addition only works on matrices of same dimension
  intersect_genes <- Reduce(intersect, lapply(seurat_list, rownames))
  seurat_list <- lapply(seurat_list, function(x) subset(x, features = intersect_genes))
  
  ## get scaled data matrices
  list_x_mats <- purrr::map(
    seurat_list,
    function(ss) {
      SeuratObject::LayerData(ss, layer = "scale.data")
    }
  )
  message(
    paste0("Processing ",
           length(list_x_mats),
           " expression matrices...")
  )
  if (is.null(names(list_x_mats))) {
    names(list_x_mats) <- paste0("X", seq_along(list_x_mats))
  }
  g_dims <- dim(list_x_mats[[1]])[1]
  l_added <- list(
    x_sharp = Matrix::Matrix(0, nrow = g_dims, ncol = g_dims, sparse = TRUE),
    x_blunt = Matrix::Matrix(0, nrow = g_dims, ncol = g_dims, sparse = TRUE),
    data_info = list()
  )
  p <- progressr::progressor(along = list_x_mats)
  for (mat_ind in seq_along(list_x_mats)) {
    x_name <- names(list_x_mats[mat_ind])
    message(paste0("Starting matrix ", x_name))
    Sys.sleep(1)
    p(message = sprintf("Processing %s", names(list_x_mats)[mat_ind]))
    x <- factorise(list_x_mats[[x_name]],
                              extrema,
                              n_threads,
                              path_to_write_angles)
    message(paste0("Adding ", x_name))
    # rows_sharp <- rownames(x[["x_sharp"]])
    # cols_sharp <- colnames(x[["x_sharp"]])
    # rows_blunt <- rownames(x[["x_blunt"]])
    # cols_blunt <- colnames(x[["x_blunt"]])
    # l_added[["x_sharp"]] <- matrixAddition(l_added[["x_sharp"]], x[["x_sharp"]])
    # l_added[["x_blunt"]] <- matrixAddition(l_added[["x_blunt"]], x[["x_blunt"]])
    l_added[["x_sharp"]] <- l_added[["x_sharp"]] + x[["x_sharp"]]
    l_added[["x_blunt"]] <- l_added[["x_blunt"]] + x[["x_blunt"]]
    l_added[["l_angles"]][[x_name]] <- x[["l_angles"]]
    l_added[["data_info"]][[x_name]] <- x[["data_info"]]
    # rownames(l_added[["x_sharp"]]) <- rows_sharp
    # colnames(l_added[["x_sharp"]]) <- cols_sharp
    # rownames(l_added[["x_blunt"]]) <- rows_blunt
    # colnames(l_added[["x_blunt"]]) <- cols_blunt
    invisible(gc())
  }

  # Implementation of the parallel processing will introduce
  # heavy memory load, which is now limited to ~4Gb of RAM for
  # an average sized dataset. Consider the tradeoffs.
  return(l_added)
}
