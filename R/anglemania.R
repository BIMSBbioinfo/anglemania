# ---------------------------------------------------------------------------
# The anglemania function extracts genes that have high biological information
# and are invariant to batch effects. It works on a SingleCellExperiment object
# and returns the same object with the selected genes and statistics
# added to the metadata.
# ---------------------------------------------------------------------------

#' @title anglemania
#' @description
#' `anglemania` computes critical angles between genes across all samples
#' provided in an \code{\link{SingleCellExperiment}}. It calculates angles,
#' transforms them to z-scores, computes statistical measures, and selects
#' genes based on specified thresholds. These genes are biologically
#' informative and invariant to batch effects.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Computes angles between genes for each batch in the
#'     \code{SingleCellExperiment} using the specified \code{method}, via
#'     \code{\link{factorise}}.
#'   \item Transforms the angles to z-scores.
#'   \item Computes statistical measures (mean z-score, signal-to-noise ratio)
#'     across batches using \code{\link{get_list_stats}}.
#'   \item Selects genes based on specified z-score thresholds using
#'     \code{\link{select_genes}}.
#' }
#'
#' The computed statistics and selected genes are added to the
#' \code{SingleCellExperiment} object, which is returned.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param batch_key Character string specifying the column in the metadata of
#' the \code{SingleCellExperiment} object that indicates which batch the cells
#' belong to.
#' @param dataset_key Character string specifying the column in the metadata of
#' the \code{SingleCellExperiment} object that indicates which dataset the cells
#' belong to. If \code{NA}, then all samples are assumed to belong to the same
#' dataset and are separated by \code{batch_key}.
#' @param min_cells_per_gene Integer specifying the minimum number of cells per
#' gene. Default is \code{1}.
#' @param min_samples_per_gene Integer specifying the minimum number of samples
#' per gene. Default is \code{2}.
#' @param allow_missing_features Logical indicating whether to allow missing
#' features. Default is \code{FALSE}.
#' @param zscore_mean_threshold Numeric value specifying the threshold
#'  for the mean z-score. Default is \code{2.5}.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   signal-to-noise z-score. Default is \code{2.5}.
#' @param max_n_genes Integer specifying the maximum number of genes to select.
#' @param method Character string specifying the method to use for calculating
#' the relationship between gene pairs. Default is \code{"cosine"}.
#' Other options include \code{"spearman"}
#' @param permute_row_or_column Character "row" or "column", whether
#' permutations
#' should be executed row-wise or column wise. Default is \code{"column"}
#' @param permutation_function Character "sample" or "permute_nonzero". If
#' sample,then sample is used for constructing background distributions. If
#' permute_nonzero, then only non-zero values are permuted. Default is
#' \code{"sample"}
#' @param prefilter_threshold Numeric value specifying the threshold
#' prefiltering
#'   genes. Speeds up gene selection.
#' @param normalization_method Character "divide_by_total_counts" or
#'   "scale_by_total_counts". Default is \code{"divide_by_total_counts"}
#' @param verbose Logical indicating whether to print progress messages.
#' @return An updated \code{SingleCellExperiment} object with computed
#'   statistics and selected genes based on the specified thresholds.
#'   The results are stored in the metadata of the \code{SingleCellExperiment}
#'   object.
#'
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @seealso
#'   \code{\link{get_list_stats}},
#'   \code{\link{select_genes}},
#'   \code{\link{factorise}},
#'   \code{\link[bigstatsr]{big_apply}},
#'   \url{https://arxiv.org/abs/1306.0256}
#'
#' @examples
#' # Set seed (optional)
#' set.seed(1)
#' sce <- sce_example()
#' sce <- anglemania(
#'   sce,
#'   batch_key = "batch",
#'   method = "cosine",
#'   zscore_mean_threshold = 2,
#'   zscore_sn_threshold = 2
#' )
#'
#' # Access the selected genes
#' selected_genes <- get_anglemania_genes(sce)
#' selected_genes[1:10]
#' @export
anglemania <- function(
    sce,
    batch_key,
    dataset_key = NA_character_,
    min_cells_per_gene = 1,
    min_samples_per_gene = 2,
    allow_missing_features = FALSE,
    zscore_mean_threshold = 2.5,
    zscore_sn_threshold = 2.5,
    max_n_genes = NULL,
    method = "cosine",
    permute_row_or_column = "column",
    permutation_function = "sample",
    prefilter_threshold = 0.5,
    normalization_method = "divide_by_total_counts",
    verbose = TRUE
) {
    # Check parameters
    S4Vectors::metadata(sce)$anglemania$params <- check_params(
        sce = sce,
        batch_key = batch_key,
        dataset_key = dataset_key,
        zscore_mean_threshold = zscore_mean_threshold,
        zscore_sn_threshold = zscore_sn_threshold,
        max_n_genes = max_n_genes,
        method = method,
        min_cells_per_gene = min_cells_per_gene,
        min_samples_per_gene = min_samples_per_gene,
        allow_missing_features = allow_missing_features,
        permute_row_or_column = permute_row_or_column,
        permutation_function = permutation_function,
        prefilter_threshold = prefilter_threshold,
        normalization_method = normalization_method,
        verbose = verbose
    )
    {
        # Process inputs
        pbapply::pboptions(
            type = "timer",
            style = 1,
            char = "=",
            title = "anglemania"
        )
    }

    if (verbose) {
        message("Preparing input...")
    }
    sce <- add_unique_batch_key(
        sce,
        dataset_key = dataset_key,
        batch_key = batch_key
    )
    S4Vectors::metadata(sce)$anglemania$weights <- .set_weights(
        col_data = SummarizedExperiment::colData(sce),
        batch_key = batch_key,
        dataset_key = dataset_key
    )
    weights <- setNames(
        S4Vectors::metadata(sce)$anglemania$weights$weight,
        S4Vectors::metadata(sce)$anglemania$weights$anglemania_batch
    )
    # split cells/barcodes by batch
    barcodes_by_batch <- split(
        rownames(SummarizedExperiment::colData(sce)),
        SummarizedExperiment::colData(sce)$anglemania_batch
    )
    # we separate the counts for each batch and convert them to FBM objects
    # later on because we will use the file-backed matrices to compute the
    #angles
    vmessage(
        verbose,
        "Filtering each batch to at least ",
        min_cells_per_gene,
        " cells per gene..."
    )
    S4Vectors::metadata(sce)$anglemania$matrix_list <-
        pbapply::pblapply(
            barcodes_by_batch,
            function(barcodes) {
                mat <- SingleCellExperiment::counts(sce)[, barcodes]
                mat <- mat[
                    Matrix::rowSums(mat > 0) >= min_cells_per_gene,
                ]
                # convert to FBM
                mat
            },
            cl = min(bigstatsr::nb_cores() - 1, 4)
        )
    names(S4Vectors::metadata(sce)$anglemania$matrix_list) <- names(
        barcodes_by_batch
    )
    S4Vectors::metadata(
        sce
    )$anglemania$intersect_genes <- get_intersect_genes(
        matrix_list = S4Vectors::metadata(sce)$anglemania$matrix_list,
        allow_missing_features = allow_missing_features,
        min_samples_per_gene = min_samples_per_gene,
        verbose = verbose
    )
    #
    S4Vectors::metadata(
        sce
    )$anglemania$matrix_list <- prepare_matrices(
        matrix_list = S4Vectors::metadata(sce)$anglemania$matrix_list,
        intersect_genes = S4Vectors::metadata(
            sce
        )$anglemania$intersect_genes,
        verbose = verbose
    )

    # compute angles and transform to z-scores
    if (verbose) {
        message("Computing angles and transforming to z-scores...")
    }
    S4Vectors::metadata(
        sce
    )$anglemania$matrix_list <- pbapply::pblapply(
        S4Vectors::metadata(sce)$anglemania$matrix_list,
        function(x) {
            factorise(
                x_mat = x,
                method = method,
                seed = 1,
                permute_row_or_column = permute_row_or_column,
                permutation_function = permutation_function,
                normalization_method = normalization_method
            )
        }
    )

    vmessage(verbose, "Computing statistics...")
    S4Vectors::metadata(sce)$anglemania$list_stats <- get_list_stats(
        matrix_list = S4Vectors::metadata(sce)$anglemania$matrix_list,
        weights = weights,
        verbose = verbose
    )
    invisible(gc())

    vmessage(verbose, "Pre-filtering features...")
    S4Vectors::metadata(
        sce
    )$anglemania$prefiltered_df <- prefilter_angl(
        snr_zscore_matrix = S4Vectors::metadata(
            sce
        )$anglemania$list_stats$sn_zscore,
        mean_zscore_matrix = S4Vectors::metadata(
            sce
        )$anglemania$list_stats$mean_zscore,
        zscore_mean_threshold = prefilter_threshold,
        zscore_sn_threshold = prefilter_threshold,
        verbose = verbose
    )

    vmessage(verbose, "Extracting filtered features...")
    sce <- select_genes(
        sce,
        zscore_mean_threshold = zscore_mean_threshold,
        zscore_sn_threshold = zscore_sn_threshold,
        max_n_genes = max_n_genes,
        verbose = verbose
    )

    return(sce)
}
