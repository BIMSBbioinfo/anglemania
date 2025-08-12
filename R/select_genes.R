# ---------------------------------------------------------------------------
# Filter and select the anglemania genes from an SCE, for which we already
# ran the anglemania function
# ---------------------------------------------------------------------------
#' @title Select genes
#' @description Select genes from a SingleCellExperiment object based on
#' mean z-score and the signal-to-noise ratio of angles between gene pairs
#' across batches.
#' @name select_genes
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
#' @describeIn select_genes Prefilter gene pairs from the mean and SNR z-scores
#' based on thresholds, to simplify downstream filtering.
#'
#' @param snr_zscore_matrix A \code{bigstatsr::FBM} object containing the SNR
#'   z-scores.
#' @param mean_zscore_matrix A \code{bigstatsr::FBM} object containing the mean
#'   z-scores.
#' @param zscore_mean_threshold Numeric value specifying the threshold for the
#'   absolute mean z-score. Default is 1.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   SNR z-score. Default is 1.
#' @return A data frame containing the prefiltered gene pairs.
#' @examples
#' library(SingleCellExperiment)
#' sce <- sce_example()
#' sce <- anglemania(sce, batch_key = "batch")
#' snr_zscore_matrix <- metadata(sce)$anglemania$list_stats$sn_zscore
#' mean_zscore_matrix <- metadata(sce)$anglemania$list_stats$mean_zscore
#' prefiltered_df <- prefilter_angl(
#'   snr_zscore_matrix,
#'   mean_zscore_matrix,
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
    snr_zscore_matrix,
    mean_zscore_matrix,
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
        BM_sn = snr_zscore_matrix,
        BM_mean = mean_zscore_matrix,
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
            BM_sn = snr_zscore_matrix,
            BM_mean = mean_zscore_matrix,
            zscore_mean_threshold = zscore_mean_threshold,
            zscore_sn_threshold = zscore_sn_threshold
        )
    }

    return(prefiltered_df)
}

# ---------------------------------------------------------------------------
#' @describeIn select_genes Select genes from the mean and SNR z-score matrices
#' stored in the SCE based on thresholds for mean and SNR.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param zscore_mean_threshold Numeric value specifying the threshold for the
#'   absolute mean z-score. Default is 2.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   SNR z-score. Default is 2.
#' @param max_n_genes Integer specifying the maximum number of genes to select.
#'   If \code{NULL}, all genes that pass the thresholds are used. Default is
#'   \code{NULL}.
#' @param direction whether to select genes with positive, negative or both
#' mean z-scores. Default is "both"
#' @param adjust_thresholds whether to automatically adjust threholds if the
#' selected genes do not meet the thresholds. Default is TRUE
#' @return The input \code{anglemaniaObject} with the
#'   \code{integration_genes} slot updated to include the selected genes and
#'   their statistical information.
#' @importFrom stats quantile
#' @importFrom dplyr filter
#' @importFrom  SummarizedExperiment rowData
#' @details
#' The function performs the following steps:
#' \enumerate{
#'  \item If \code{max_n_genes} is not specified, it uses all genes that pass
#'     the thresholds.
#'  \item Identifies gene pairs where both the mean z-score and SNR z-score
#'     exceed the specified thresholds.
#'  \item If no gene pairs meet the criteria, it adjusts the thresholds to the
#'     99th percentile values of the corresponding statistics and re-selects.
#'  \item Extracts unique genes from the selected gene pairs using
#'     \code{\link{extract_rows_for_unique_genes}}.
#' }
#' @examples
#' sce <- sce_example()
#' sce <- anglemania(
#'   sce,
#'   batch_key = "batch",
#'   zscore_mean_threshold = 2.5,
#'   zscore_sn_threshold = 2.5
#' )
#' anglemania_genes <- get_anglemania_genes(sce)
#' # View the selected genes and use for integration
#' head(anglemania_genes)
#' length(anglemania_genes)
#' # Adjust thresholds to 2 and select genes
#' sce <- select_genes(
#'   sce,
#'   zscore_mean_threshold = 2,
#'   zscore_sn_threshold = 2
#' )
#' anglemania_genes <- get_anglemania_genes(sce)
#' head(anglemania_genes)
#' length(anglemania_genes)
#' @seealso \code{\link{extract_rows_for_unique_genes}},
#'   \code{\link{get_intersect_genes}}, \code{\link{get_list_stats}}
#' @export
select_genes <- function(
    sce,
    zscore_mean_threshold = 2,
    zscore_sn_threshold = 2,
    max_n_genes = NULL,
    direction = "both",
    adjust_thresholds = TRUE,
    verbose = TRUE
) {
    if (!direction %in% c("both", "positive", "negative")) {
        stop("direction can only be 'both', 'positive' or 'negative'")
    }
    # check if anglemania element is in metadata
    if (!"anglemania" %in% names(S4Vectors::metadata(sce))) {
        stop("please run anglemania first")
    }
    prefiltered_df <- S4Vectors::metadata(
        sce
    )$anglemania$prefiltered_df
    prefiltered_df$geneA <-
        S4Vectors::metadata(sce)$anglemania$intersect_genes[
            prefiltered_df$geneA
        ]
    prefiltered_df$geneB <-
        S4Vectors::metadata(sce)$anglemania$intersect_genes[
            prefiltered_df$geneB
        ]
    if (is.null(max_n_genes)) {
        # If no max_n_genes specified, use all genes that pass the threshold
        max_n_genes <- length(
            S4Vectors::metadata(sce)$anglemania$intersect_genes
        )
    }

    # Selects the direction of conserved genes
    if (direction == "both") {
        filtered_genes_df <- subset(
            prefiltered_df,
            sn_zscore >= zscore_sn_threshold &
                abs(prefiltered_df$mean_zscore) >=
                    zscore_mean_threshold
        )
    } else if (direction == "positive") {
        filtered_genes_df <- subset(
            prefiltered_df,
            sn_zscore >= zscore_sn_threshold &
                prefiltered_df$mean_zscore >= zscore_mean_threshold
        )
    } else if (direction == "negative") {
        filtered_genes_df <- subset(
            prefiltered_df,
            sn_zscore >= zscore_sn_threshold &
                prefiltered_df$mean_zscore <= zscore_mean_threshold
        )
    }

    # Adjust thresholds if no genes passed the cutoff
    while (nrow(filtered_genes_df) == 0 && adjust_thresholds) {
        vmessage(verbose, "No genes passed the cutoff.")
        zscore_sn_threshold <- zscore_sn_threshold - 0.1
        zscore_mean_threshold <- zscore_mean_threshold - 0.1
        vmessage(
            verbose,
            paste0(
                "Decreasing zscore mean and sn thresholds by 0.1: \n",
                "zscore_mean_threshold: ",
                round(zscore_mean_threshold, 2),
                "\n",
                "zscore_sn_threshold: ",
                round(zscore_sn_threshold, 2)
            )
        )
        filtered_genes_df <- prefiltered_df |>
            dplyr::filter(
                abs(mean_zscore) >= zscore_mean_threshold &
                    sn_zscore >= zscore_sn_threshold
            )
    }

    # Order data frame
    filtered_genes_df <- filtered_genes_df[
        order(abs(filtered_genes_df$mean_zscore), decreasing = TRUE),
    ]
    # Extract the unique genes from the gene pairs
    # that passed the thresholds
    anglemania_genes <- extract_rows_for_unique_genes(
        filtered_genes_df,
        max_n_genes
    )

    # Assign the anglemania genes to the metadata of SCE/SE object
    # and the rowData of SCE/SE object
    S4Vectors::metadata(
        sce
    )$anglemania$anglemania_genes <- anglemania_genes
    SummarizedExperiment::rowData(sce)$anglemania_genes <-
        rownames(sce) %in% anglemania_genes

    # update anglemania params slot with the arguments
    S4Vectors::metadata(
        sce
    )$anglemania$params$zscore_mean_threshold <-
        zscore_mean_threshold
    S4Vectors::metadata(sce)$anglemania$params$zscore_sn_threshold <-
        zscore_sn_threshold
    S4Vectors::metadata(sce)$anglemania$params$direction <- direction
    S4Vectors::metadata(
        sce
    )$anglemania$params$max_n_genes <- max_n_genes
    vmessage(
        verbose,
        "Selected ",
        length(anglemania_genes),
        " genes for integration."
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
