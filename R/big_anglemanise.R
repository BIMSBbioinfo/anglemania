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
#' @importFrom bigstatsr big_transpose
#' @importFrom bigstatsr big_apply
#' @importFrom bigstatsr as_FBM
#' @importFrom margittr %>%
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
#' @param seurat_list seurat list of seurat objects with scaled data.
#' @param extrema double. Fraction of the angles
#'   to be cut from both sides of an approximated angle
#'   distribution.
#' @param n_cores integer. Number of cores to use for
#'   the computation of the cosine distances.
#' @return list. First two elements are sparse matrices
#'   recording the number of sharp and blunt angles across
#'   the datasets. Third and fourth elements are lists with
#'   angles statistics and paths to the values of critical
#'   angles.
#' @seealso https://arxiv.org/abs/1306.0256
#' @export big_anglemanise
big_anglemanise <- function(seurat_list, # nolint
                            extrema = 0.005,
                            n_cores = 4) {
    ############## Validate inputs ###########################

    if (!is.list(seurat_list) || any(sapply(seurat_list, function(x) attr(class(x), "package") != "SeuratObject"))) {
        stop("seurat_list needs to be a list of Seurat objects")
    }
    if (!is.numeric(n_cores) || n_cores < 1) {
        stop("n_cores has to be a positive integer")
    }
    if (length(extrema) > 2) {
        stop("extrema should be numeric of length either 1 or 2")
    }
    if (length(extrema) == 1) {
        extrema <- c(extrema, 1 - extrema)
    }


    ############## Process inputs ###########################
    # Get the list of count matrices
    pboptions(
        type = "timer",
        style = 1,
        char = "=",
        title = "extract count matrices"
    )
    message("Extracting count matrices...")
    list_x_mats <- pbapply::pblapply(seurat_list, function(x) {
        x <- SeuratObject::GetAssayData(x, assay = "RNA")
    }, cl = n_cores)

    # Assign names to the list of count matrices
    # If names(seurat_list) is NULL, it assigns sequential names
    # based on the position of matrices in the list.
    # Otherwise, it assigns basename of the corresponding
    # seurat object names.
    if (is.null(names(list_x_mats))) {
        names(list_x_mats) <- paste0("X", seq_along(list_x_mats))
    } else {
        names(list_x_mats) <- sapply(names(seurat_list), basename)
    }


    # Calculate correlation matrix of genes for each sample and only keep significant ones
    #   for details see big_factorise
    # l_processed <- parallel::mclapply(names(list_x_mats), function(name) {
    #     x <- list_x_mats[[name]]
    #     x <- big_factorise(
    #         x_mat = x,
    #         name = name,
    #         extrema
    #     )
    # }, mc.cores = n_cores)


    pbapply::pboptions(
        type = "timer",
        style = 1,
        char = "=",
        title = "big_anglemanise"
    )
    message("Calculating significant correlations...")
    l_processed <- pbapply::pblapply(names(list_x_mats), function(name) {
        x <- list_x_mats[[name]]
        x <- big_factorise(
            x_mat = x,
            name = name,
            extrema
        )
    }, cl = n_cores)

    gc()
    return(l_processed)
}
