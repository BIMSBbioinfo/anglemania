# ---------------------------------------------------------------------------
# The anglemania function extracts genes that have high biological information
# and are invariant to batch effects. It works on a SingleCellExperiment object
# and returns the same object with the selected genes and statistics
# added to the metadata.
# ---------------------------------------------------------------------------

#' @title anglemania
#' @description
#' `anglemania` computes critical angles between genes across all samples
#' provided in an
#' \link[SingleCellExperiment:SingleCellExperiment-class]{SingleCellExperiment}
#' object. It calculates angles, transforms them to z-scores, computes
#' statistical measures, and selects the top genes based on mean and standard
#' deviation of z-scores.
#' These genes are biologically informative and invariant to batch effects.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Optionally, when \code{use_simpsons = TRUE}, groups cells into
#'     microclusters within each batch and derives per-cell weights so that
#'     each microcluster contributes equally to the correlation estimate.
#'   \item Computes angles between genes for each batch in the
#'     \code{SingleCellExperiment} using the specified \code{method}, via
#'     \code{\link{factorise}}.
#'   \item Transforms the angles to z-scores.
#'   \item Computes statistical measures (mean z-score, signal-to-noise ratio)
#'     across batches using \code{\link{get_list_stats}}.
#'   \item Selects the top n genes based on mean and standard deviation of
#'     z-scores using \code{\link{select_genes}}.
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
#' @param max_n_genes Integer specifying the maximum number of genes to select.
#' @param min_cells_per_gene Integer specifying the minimum number of cells per
#' gene. Default is \code{1}.
#' @param min_samples_per_gene Integer specifying the minimum number of samples
#' per gene. Default is \code{2}.
#' @param allow_missing_features Logical indicating whether to allow missing
#' features. Default is \code{FALSE}.
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
#' @param use_simpsons Logical indicating whether to apply Simpson's paradox
#' correction. When \code{TRUE}, cells are grouped into microclusters within
#' each batch and gene-gene correlations are computed with per-cell weights
#' \code{1 / |cluster|}, so that each microcluster contributes equally
#' regardless of size. This removes correlation signal that is driven purely by
#' differing cell type proportions across batches. Default is \code{FALSE}.
#' @param simpsons_n_clusters Integer nominally specifying the target number
#' of microclusters per batch. \strong{This argument currently has no effect
#' for values of 20 or below.} It sets the Leiden resolution to
#' \code{max(1, simpsons_n_clusters / 20)}, so every value from 2 to 20 maps to
#' resolution 1 and yields identical clustering; above 20 the resolution changes
#' but \code{igraph::cluster_leiden(objective_function = "modularity")}
#' responds to it only weakly and non-monotonically. In practice the clustering
#' resolves whole cell types rather than finer microclusters. The published
#' benchmarks were produced at the default and are unaffected, but do not
#' expect to tune this. Default is \code{15}.
#' @param simpsons_n_pcs Integer specifying the number of principal components
#' used for the microclustering embedding. Default is \code{20}.
#' @param simpsons_min_cells Integer specifying the minimum number of cells a
#' microcluster may contain. Clusters below this size are iteratively merged
#' into their nearest neighbour. Default is \code{30}.
#' @param simpsons_min_cluster_size Integer specifying the minimum cluster size
#' for a cell to receive a non-zero weight. Cells in smaller clusters are
#' excluded from the correlation estimate. Default is \code{20}.
#' @param simpsons_min_cells_cor Integer specifying the minimum number of cells
#' a microcluster must contain to contribute a within-cluster correlation
#' during the post-hoc analysis. Only used when
#' \code{simpsons_weight_in_ranking > 0}. Default is \code{20}.
#' @param simpsons_top_quantile Numeric in \code{[0, 1]} specifying the
#' quantile of within-cluster correlations used to summarise a gene pair's
#' peak local correlation. Only used when
#' \code{simpsons_weight_in_ranking > 0}. Default is \code{0.95}.
#' @param simpsons_weight_in_ranking Numeric in \code{[0, 1]} specifying how
#' much weight the Simpson's preservation score receives as a third criterion
#' in \code{\link{select_genes}}, alongside mean z-score and z-score standard
#' deviation. The two base criteria are rescaled to sum to
#' \code{1 - simpsons_weight_in_ranking}. Setting this above \code{0} triggers
#' the post-hoc local correlation analysis, which is computationally
#' expensive. Default is \code{0} (disabled).
#' @param simpsons_use_binary_pca Logical indicating whether to build the
#' clustering embedding from the binarised detection matrix rather than from
#' normalized expression. Benchmarks show this costs 0.4-0.8\% DE precision;
#' retained as an option but not recommended. Default is \code{FALSE}.
#' @param simpsons_use_density_weights Logical indicating whether to replace
#' the discrete cluster-inverse weights with continuous kNN density-based
#' weights. Benchmarks show this costs 3-7\% DE precision; retained as an
#' option but not recommended. Default is \code{FALSE}.
#' @param simpsons_density_k Integer specifying the number of nearest
#' neighbours for density estimation. Only used when
#' \code{simpsons_use_density_weights = TRUE}. Default is \code{30}.
#' @param simpsons_use_count_split Logical indicating whether to use count
#' splitting (Neufeld et al., \url{https://arxiv.org/abs/2307.12985}) so that
#' clustering and correlation estimation use independent Poisson folds of the
#' counts. This avoids double-dipping and is recommended on real data, where
#' uncorrected Simpson's weighting tends to overcorrect. Requires the
#' \pkg{countsplit} package. Default is \code{FALSE}.
#' @param simpsons_count_split_prop Numeric in \code{[0.1, 0.9]} specifying the
#' proportion of counts allocated to the first fold (used for PCA and
#' clustering); the remainder is used for correlations. Only used when
#' \code{simpsons_use_count_split = TRUE}. Default is \code{0.5}.
#' @return An updated \code{SingleCellExperiment} object with computed
#'   statistics and selected genes.
#'   The results are stored in the metadata of the \code{SingleCellExperiment}
#'   object.
#'
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
#' @importFrom stats setNames
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
#'   method = "cosine"
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
    max_n_genes = 2000,
    min_cells_per_gene = 1,
    min_samples_per_gene = 2,
    allow_missing_features = FALSE,
    method = "cosine",
    permute_row_or_column = "column",
    permutation_function = "sample",
    prefilter_threshold = 0.5,
    normalization_method = "divide_by_total_counts",
    verbose = TRUE,
    use_simpsons = FALSE,
    simpsons_n_clusters = 15,
    simpsons_n_pcs = 20,
    simpsons_min_cells = 30,
    simpsons_min_cluster_size = 20,
    simpsons_min_cells_cor = 20,
    simpsons_top_quantile = 0.95,
    simpsons_weight_in_ranking = 0.0,
    simpsons_use_binary_pca = FALSE,
    simpsons_use_density_weights = FALSE,
    simpsons_density_k = 30,
    simpsons_use_count_split = FALSE,
    simpsons_count_split_prop = 0.5
) {
    # Check parameters
    S4Vectors::metadata(sce)$anglemania$params <- check_params(
        sce = sce,
        batch_key = batch_key,
        dataset_key = dataset_key,
        max_n_genes = max_n_genes,
        method = method,
        min_cells_per_gene = min_cells_per_gene,
        min_samples_per_gene = min_samples_per_gene,
        allow_missing_features = allow_missing_features,
        permute_row_or_column = permute_row_or_column,
        permutation_function = permutation_function,
        prefilter_threshold = prefilter_threshold,
        normalization_method = normalization_method,
        verbose = verbose,
        use_simpsons = use_simpsons,
        simpsons_n_clusters = simpsons_n_clusters,
        simpsons_n_pcs = simpsons_n_pcs,
        simpsons_min_cells = simpsons_min_cells,
        simpsons_min_cluster_size = simpsons_min_cluster_size,
        simpsons_min_cells_cor = simpsons_min_cells_cor,
        simpsons_top_quantile = simpsons_top_quantile,
        simpsons_weight_in_ranking = simpsons_weight_in_ranking,
        simpsons_use_binary_pca = simpsons_use_binary_pca,
        simpsons_use_density_weights = simpsons_use_density_weights,
        simpsons_density_k = simpsons_density_k,
        simpsons_use_count_split = simpsons_use_count_split,
        simpsons_count_split_prop = simpsons_count_split_prop
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
    S4Vectors::metadata(
        sce
    )$anglemania$params$dataset_weights <- .set_weights(
        col_data = SummarizedExperiment::colData(sce),
        batch_key = batch_key,
        dataset_key = dataset_key
    )
    params <- S4Vectors::metadata(sce)$anglemania$params
    dataset_weights <- setNames(
        params$dataset_weights$weight,
        params$dataset_weights$anglemania_batch
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
    # Simpson's paradox correction: microcluster before FBM conversion
    simpsons_data <- NULL
    original_matrices <- NULL
    if (use_simpsons) {
        vmessage(verbose, "Preparing Simpson's correction...")
        simpsons_data <- prepare_simpsons(
            matrix_list = S4Vectors::metadata(sce)$anglemania$matrix_list,
            n_clusters = simpsons_n_clusters,
            n_pcs = simpsons_n_pcs,
            min_cells = simpsons_min_cells,
            min_cluster_size = simpsons_min_cluster_size,
            use_binary_pca = simpsons_use_binary_pca,
            use_density_weights = simpsons_use_density_weights,
            density_k = simpsons_density_k,
            use_count_split = simpsons_use_count_split,
            count_split_prop = simpsons_count_split_prop,
            verbose = verbose
        )
        S4Vectors::metadata(sce)$anglemania$simpsons_data <- simpsons_data

        # If count splitting, replace correlation matrices with fold 2
        if (simpsons_use_count_split &&
            !is.null(simpsons_data$cor_matrices)) {
            S4Vectors::metadata(sce)$anglemania$matrix_list <-
                simpsons_data$cor_matrices
        }
        # Keep original sparse matrices for posthoc local correlations
        if (simpsons_weight_in_ranking > 0) {
            original_matrices <- S4Vectors::metadata(
                sce
            )$anglemania$matrix_list
        }
    }

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
    fbm_list <- S4Vectors::metadata(sce)$anglemania$matrix_list
    batch_names <- names(fbm_list)
    S4Vectors::metadata(
        sce
    )$anglemania$matrix_list <- pbapply::pblapply(
        stats::setNames(batch_names, batch_names),
        function(bn) {
            w <- if (use_simpsons) {
                simpsons_data$weights_per_batch[[bn]]
            } else {
                NULL
            }
            factorise(
                x_mat = fbm_list[[bn]],
                method = method,
                seed = 1,
                permute_row_or_column = permute_row_or_column,
                permutation_function = permutation_function,
                normalization_method = normalization_method,
                cell_weights = w
            )
        }
    )

    vmessage(verbose, "Computing statistics...")
    S4Vectors::metadata(sce)$anglemania$list_stats <- get_list_stats(
        matrix_list = S4Vectors::metadata(sce)$anglemania$matrix_list,
        weights = dataset_weights,
        verbose = verbose
    )
    invisible(gc())

    vmessage(verbose, "Pre-filtering features...")
    sce <- prefilter_angl(
        sce,
        zscore_mean_threshold = prefilter_threshold,
        zscore_sn_threshold = prefilter_threshold,
        verbose = verbose
    )

    # Post-hoc Simpson's local correlation analysis
    if (use_simpsons && simpsons_weight_in_ranking > 0) {
        vmessage(verbose, "Computing Simpson's local correlations...")
        sce <- posthoc_simpsons(
            sce = sce,
            original_matrices = original_matrices,
            min_cells_cor = simpsons_min_cells_cor,
            top_quantile = simpsons_top_quantile,
            verbose = verbose
        )
    }

    vmessage(verbose, "Extracting filtered features...")
    sce <- select_genes(
        sce,
        max_n_genes = max_n_genes,
        verbose = verbose
    )

    return(sce)
}
