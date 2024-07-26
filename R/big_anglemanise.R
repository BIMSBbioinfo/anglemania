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
#' @importFrom magrittr %>%
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
#' @param anglem_object seurat object. The Seurat object should be a combined seurat object of all the datasets and samples.
#' @param fdr_threshold double. The FDR threshold to apply
#'   to the q-values. Defaults to 0.001.
#' @param n_cores integer. Number of cores to use for
#'   the computation of the cosine distances.
#' @return list. First two elements are sparse matrices
#'   recording the number of sharp and blunt angles across
#'   the datasets. Third and fourth elements are lists with
#'   angles statistics and paths to the values of critical
#'   angles.
#' @seealso https://arxiv.org/abs/1306.0256
#' @export big_anglemanise
big_anglemanise <- function(anglem_object, # nolint
                            zscore_mean_threshold = 2,
                            zscore_sn_threshold  = 2,
                            max_n_genes = 2000,
                            n_cores = 3) {
    ############## Validate inputs ###########################

    if (class(anglem_object) != "anglem") {
        stop("anglem_object needs to be a anglem object")
    }

    # dataset_key and batch_key are checked in create_anglem!

    if (!is.numeric(zscore_mean_threshold) || zscore_mean_threshold < 0) {
        stop("zscore_mean_threshold has to be a non-negative number")
    }
    if (!is.numeric(zscore_sn_threshold) || zscore_sn_threshold < 0) {
        stop("zscore_sn_threshold has to be a non-negative number")
    }
    if (!is.numeric(max_n_genes) || max_n_genes < 1) {
        stop("max_n_genes has to be a positive integer")
    }
    if (!is.numeric(n_cores) || n_cores < 1) {
        stop("n_cores has to be a positive integer")
    }
    ############## Process inputs ###########################
    pbapply::pboptions(
        type = "timer",
        style = 1,
        char = "=",
        title = "big_anglemanise"
    )
    message("Computing correlations and transforming to z-scores...")
    matrix_list(anglem_object) <- pbapply::pblapply(matrix_list(anglem_object), function(x) {
        x <- big_factorise(x_mat = x,
                           seed = 1)
    })
    
    message("Computing statistics...")
    list_stats(anglem_object) <- get_list_stats(anglem_object)
    invisible(gc())

    message("Filtering features...")
    anglem_object <- select_genes(
        anglem_object,
        zscore_mean_threshold = zscore_mean_threshold,
        zscore_sn_threshold = zscore_sn_threshold,
        max_n_genes = max_n_genes
    )
    
    return(anglem_object)
}
