# ---------------------------------------------------------------------------
# Filter and select the anglemania genes from an SCE, for which we already
# ran the anglemania function
# ---------------------------------------------------------------------------
#' @title Select genes from an anglemania-processed SCE
#' @description Select genes from a SingleCellExperiment object based on
#' mean z-score and the signal-to-noise ratio of angles between gene pairs
#' across batches.
#' @name select_genes
#' @rdname select_genes
#' @keywords internal
NULL


#' @describeIn select_genes Prefilter gene pairs from the mean and SNR z-scores
#' based on thresholds, to simplify downstream filtering.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param zscore_mean_threshold Numeric value specifying the threshold for the
#'   absolute mean z-score. Default is 1.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   SNR z-score. Default is 1.
#' @param verbose Logical value indicating whether to print progress messages.
#'   Default is \code{TRUE}.
#' @return A data frame containing the prefiltered gene pairs.
#' @examples
#' library(SingleCellExperiment)
#' sce <- sce_example()
#' sce <- anglemania(sce, batch_key = "batch")
#' prefiltered_df <- prefilter_angl(
#'   sce,
#'   zscore_mean_threshold = 1,
#'   zscore_sn_threshold = 1
#' )
#' head(prefiltered_df)
#' @details
#' The function performs the following steps:
#' \enumerate{
#'  \item Identifies gene pairs where both the mean z-score and SNR z-score
#'     exceed the specified thresholds.
#' }
#' @useDynLib anglemania, .registration = TRUE
#' @export
prefilter_angl <- function(
    sce,
    zscore_mean_threshold = 1,
    zscore_sn_threshold = 1,
    verbose = TRUE
) {
    if (zscore_mean_threshold <= 0 || zscore_sn_threshold <= 0) {
        stop(
            "zscore_mean_threshold and zscore_sn_threshold need to be positive"
        )
    }
    prefiltered_df <- select_genes_cpp(
        BM_sn = S4Vectors::metadata(sce)$anglemania$list_stats$sn_zscore,
        BM_mean = S4Vectors::metadata(sce)$anglemania$list_stats$mean_zscore,
        BM_sd = S4Vectors::metadata(sce)$anglemania$list_stats$sds_zscore,
        zscore_mean_threshold = zscore_mean_threshold,
        zscore_sn_threshold = zscore_sn_threshold
    )
    while (nrow(prefiltered_df) == 0) {
        if (zscore_mean_threshold <= 0 || zscore_sn_threshold <= 0) {
            stop(
                "zscore_mean_threshold and zscore_sn_threshold ",
                "need to be positive"
            )
        }
        vmessage(
            verbose,
            "No genes passed the cutoff. Decreasing thresholds by 0.1..."
        )
        zscore_mean_threshold <- zscore_mean_threshold - 0.1
        zscore_sn_threshold <- zscore_sn_threshold - 0.1
        prefiltered_df <- select_genes_cpp(
            BM_sn = S4Vectors::metadata(sce)$anglemania$list_stats$sn_zscore,
            BM_mean = S4Vectors::metadata(sce)$anglemania$list_stats$mean_zscore,
            BM_sd = S4Vectors::metadata(sce)$anglemania$list_stats$sds_zscore,
            zscore_mean_threshold = zscore_mean_threshold,
            zscore_sn_threshold = zscore_sn_threshold
        )
    }
    prefiltered_df$geneA <-
        S4Vectors::metadata(sce)$anglemania$intersect_genes[
            prefiltered_df$geneA
        ]
    prefiltered_df$geneB <-
        S4Vectors::metadata(sce)$anglemania$intersect_genes[
            prefiltered_df$geneB
        ]
    S4Vectors::metadata(sce)$anglemania$prefiltered_df <- prefiltered_df
    return(sce)
}

# ---------------------------------------------------------------------------
#' @describeIn select_genes Select the top n genes on the weighted sum
#' of the ranks of the mean z-score and SNR z-score of the gene pairs.
#' 
#' @param sce A \code{SingleCellExperiment} object.
#' @param max_n_genes Integer specifying the maximum number of genes to select.
#' If \code{NULL}, all genes that passed the prefiltering thresholds are used.
#' Default is \code{2000}.
#' @param score_weights A vector of two numeric values specifying the weights
#' for the mean z-score and standard deviation of z-score, respectively.
#' Default is \code{c(0.4, 0.6)} for a greater emphasis on the standard
#' deviation of z-score.
#' @return The input \code{SingleCellExperiment} object with the
#'   \code{anglemania_genes} slot updated to include the selected genes and
#'   their statistical information.
#' @importFrom dplyr filter
#' @importFrom  SummarizedExperiment rowData
#' @details
#' Selects the top n genes based on the weighted sum of the ranked mean
#' and standard deviation of the z-score of the correlations between gene pairs.
#' @examples
#' sce <- sce_example()
#' sce <- anglemania(
#'     sce,
#'     batch_key = "batch",
#'     max_n_genes = 20
#' )
#' anglemania_genes <- get_anglemania_genes(sce)
#' # View the selected genes and use for integration
#' head(anglemania_genes)
#' length(anglemania_genes)
#' sce <- select_genes(
#'     sce,
#'     max_n_genes = 10
#' )
#' anglemania_genes <- get_anglemania_genes(sce)
#' head(anglemania_genes)
#' length(anglemania_genes)
#' @seealso \code{\link{extract_rows_for_unique_genes}},
#'   \code{\link{get_intersect_genes}}, \code{\link{get_list_stats}}
#' @export
select_genes <- function(
    sce,
    max_n_genes = 2000,
    score_weights = c(0.4, 0.6),
    verbose = TRUE
) {
    if (!"anglemania" %in% names(S4Vectors::metadata(sce))) {
        stop("please run anglemania first")
    }
    if (!checkmate::test_integerish(
        max_n_genes,
        upper = length(S4Vectors::metadata(sce)$anglemania$intersect_genes)
    )) {
        vmessage(
            verbose,
            max_n_genes, " is larger than the number of intersected genes.",
            "Setting max_n_genes to ",
            length(S4Vectors::metadata(sce)$anglemania$intersect_genes)
        )
        max_n_genes <-
            length(S4Vectors::metadata(sce)$anglemania$intersect_genes)
    }
    prefiltered_df <- S4Vectors::metadata(sce)$anglemania$prefiltered_df
# Check if score_weights is NULL or valid
    check_weights <- if (is.null(score_weights)) {
        TRUE
    } else {
        checkmate::check_numeric(
            score_weights,
            lower = 0,
            upper = 1,
            len = 2,
            null.ok = FALSE
        )
    }
    if (isTRUE(check_weights)) {
        prefiltered_df$rank <- prefiltered_df |>
            dplyr::mutate(
                rank_mean_zscore = rank(
                    -abs(mean_zscore),
                    ties.method = "min"
                ),
                rank_sd_zscore = rank(sd_zscore, ties.method = "min"),
                rank = rank(
                    rank_mean_zscore * score_weights[1] +
                        rank_sd_zscore * score_weights[2]
                )
            ) |>
            dplyr::pull(rank)
    } else {
        stop(check_weights) # if weights are set incorrectly, it prints
        # the error message
    }
    
    # Extract the unique genes from the gene pairs
    anglemania_genes <- extract_rows_for_unique_genes(
        prefiltered_df |> dplyr::arrange(rank),
        max_n_genes
    )

    # add metadata to SCE/SE object
    S4Vectors::metadata(sce)$anglemania$anglemania_genes <-
        anglemania_genes
    SummarizedExperiment::rowData(sce)$anglemania_genes <-
        rownames(sce) %in% anglemania_genes
    S4Vectors::metadata(sce)$anglemania$prefiltered_df <- prefiltered_df

    # update anglemania params slot
    S4Vectors::metadata(sce)$anglemania$params$max_n_genes <- max_n_genes
    S4Vectors::metadata(sce)$anglemania$params$score_weights <- score_weights
    vmessage(
        verbose,
        "Selected ", length(anglemania_genes), " genes for integration."
    )
    return(sce)
}

# ---------------------------------------------------------------------------
#' @describeIn select_genes Extract unique gene identifiers
#' from gene pairs, returning up to a specified maximum number.
#'
#' @param dt A data frame containing gene pairs, with columns \code{geneA}
#'   and \code{geneB}.
#' @param max_n_genes An integer specifying the maximum number of unique genes
#'   to return.
#' @return A vector of unique gene identifiers.
#' @details
#' The function combines the \code{geneA} and \code{geneB} columns, extracts
#' unique gene names, and returns the first \code{max_n_genes} genes. If
#' \code{max_n_genes} exceeds the number of unique genes available, all unique
#' genes are returned.
#' @examples
#' gene_pairs <- data.frame(
#'   geneA = c("Gene1", "Gene2", "Gene3", "Gene4"),
#'   geneB = c("Gene3", "Gene4", "Gene5", "Gene6")
#' )
#' unique_genes <- extract_rows_for_unique_genes(
#'   gene_pairs,
#'   max_n_genes = 3
#' )
#' print(unique_genes)
#' @seealso \code{\link{select_genes}}
#' @export
extract_rows_for_unique_genes <- function(dt, max_n_genes) {
    unique_genes <- unique(as.vector(rbind(dt$geneA, dt$geneB)))
    max_genes <- min(max_n_genes, length(unique_genes))
    unique_genes <- unique_genes[seq_len(max_genes)]
    return(unique_genes)
}
