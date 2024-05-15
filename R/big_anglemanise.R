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
#' @param seurat_list seurat list of seurat objects with scaled data.
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
big_anglemanise <- function(seurat_list, # nolint
                            zscore_mean_threshold = 2,
                            zscore_sn_threshold  = 2,
                            max_n_genes = 2000,
                            n_cores = 4) {
    ############## Validate inputs ###########################

    if (!is.list(seurat_list) || any(sapply(seurat_list, function(x) attr(class(x), "package") != "SeuratObject"))) {
        stop("seurat_list needs to be a list of Seurat objects")
    }
    if (!is.numeric(n_cores) || n_cores < 1) {
        stop("n_cores has to be a positive integer")
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
    intersect_genes <- Reduce(intersect, lapply(seurat_list, rownames))
    message("Using intersection of genes. Number of genes: ", length(intersect_genes))
    list_x_mats <- pbapply::pblapply(seurat_list, function(x) {
        x <- subset(x, features = intersect_genes)
        x <- SeuratObject::LayerData(x, layer = "counts", assay = "RNA") # GetAssayData or LayerData from SeuratObject?
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

    pbapply::pboptions(
        type = "timer",
        style = 1,
        char = "=",
        title = "big_anglemanise"
    )
    message("Computing correlations and transforming to z-scores...")
    list_zscore_mats <- pbapply::pblapply(names(list_x_mats), function(name) {
        x <- list_x_mats[[name]]
        x <- big_factorise(x_mat = x, name = name)
    }, cl = n_cores)
    
    message("Computing statistics...")
    list_stats <- get_list_stats(list_zscore_mats)
    invisible(gc())

    l_out <- list(
        list_stats = list_stats,
        gene_names = NULL
    )

    message("Filtering features...")
    l_out <- select_genes(
        l_out,
        zscore_mean_threshold = 2,
        zscore_sn_threshold = zscore_sn_threshold,
        max_n_genes = max_n_genes,
        intersect_genes = intersect_genes
    )
    
    return(l_out)
}
