# ---------------------------------------------------------------------------
# Utility functions for the anglemania package
# ---------------------------------------------------------------------------
#' @title Utility Functions for the anglemania Package
#' @description A collection of utility functions used within the
#' \pkg{anglemania} package
#' for manipulating FBMs, calculating statistics, and selecting genes.
#' @name anglemania_utils
#' @rdname anglemania_utils
#' @keywords internal
NULL

#' @describeIn anglemania_utils Convert a sparse matrix into
#' a file-backed matrix
#' (\code{\link[bigstatsr]{FBM}}) with efficient memory usage.
#'
#' @param s_mat A sparse matrix.
#' @return An \code{\link[bigstatsr]{FBM}} object from the \pkg{bigstatsr}
#'   package.
#' @importFrom bigstatsr FBM
#' @examples
#' s_mat <- Matrix::rsparsematrix(nrow = 10, ncol = 5, density = 0.3)
#' fbm_mat <- sparse_to_fbm(s_mat)
#' fbm_mat
#' @export
sparse_to_fbm <- function(s_mat) {
    n <- nrow(s_mat)
    p <- ncol(s_mat)
    X <- bigstatsr::FBM(n, p)

    bigstatsr::big_apply(
        X,
        a.FUN = function(X, ind) {
            X[, ind] <- s_mat[, ind] |> as.matrix()
            NULL
        },
        a.combine = "c",
        block.size = 1000
    )

    return(X)
}


# ---------------------------------------------------------------------------
#' @describeIn anglemania_utils replace Nan and Inf values with NA
#'
#' @description
#' Replaces all NaN and Inf values in a numeric vector with NA.
#' @param v A numeric vector.
#' @return A numeric vector with NaN and Inf values replaced with NA.
#' @examples
#' v <- c(1, 2, 3, 4, 5, 6, 7, Inf, 9, NA)
#' v <- replace_with_na(v)
#' v
#' @export
replace_with_na <- function(v) {
    v[is.nan(v)] <- NA
    v[is.infinite(v)] <- NA
    v
}


# ---------------------------------------------------------------------------
#' @describeIn anglemania_utils normalize matrix
#' Normalize a Filebacked Big Matrix (`FBM`) using either total-count
#' scaling or residuals from a linear model. Intended for use with
#' single-cell RNA-seq gene expression data.
#'
#' @param x_mat A `bigstatsr::FBM` object containing the matrix to normalize
#'   (typically genes x cells).
#' @param normalization_method A character string specifying the normalization
#'   method to use. One of `"divide_by_total_counts"` (default) or
#'   `"find_residuals"`.
#'   - `"divide_by_total_counts"` normalizes each cell by its total expression
#'     count and applies log1p.
#'   - `"find_residuals"` computes log1p-transformed residuals after regressing
#'     out total expression.
#'
#' @return The input `FBM` object with normalized values written back in place.
#'   This function modifies the input `x_mat` by reference.
#'
#' @importFrom bigstatsr big_apply
#' @importFrom checkmate assertChoice assertClass
#' @examples
#' library(bigstatsr)
#' set.seed(42)
#' mat <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
#' fbm <- as_FBM(mat)
#'
#' normalize_matrix(fbm, normalization_method = "divide_by_total_counts")[1:5, 1:5]
#' normalize_matrix(fbm, normalization_method = "find_residuals")[1:5, 1:5]
#'
#' @export
normalize_matrix <- function(
    x_mat,
    normalization_method = "divide_by_total_counts",
    verbose = TRUE
) {
    checkmate::assertChoice(
        normalization_method,
        c("divide_by_total_counts", "find_residuals")
    )
    checkmate::assertClass(x_mat, "FBM")

    bigstatsr::big_apply(x_mat, a.FUN = function(X, ind) {
        X.sub <- X[, ind, drop = FALSE]
        if (normalization_method == "divide_by_total_counts") {
            # Normalize the data:
            #   divide gene counts by the number of total counts per cell,
            #   multiply by 10,000 (scaling factor like in Seurat)
            # +1 prevents NaN values when working with partial features

            X.sub <- t(t(X.sub) / (colSums(X.sub)) * 10000)
            X.res <- log1p(X.sub)
        } else if (normalization_method == "find_residuals") {
            total <- log1p(colSums(X.sub))
            X.sub <- log1p(X.sub)
            X.res <- t(apply(X.sub, 1, function(x) residuals(lm(x ~ total))))
        }

        X[, ind] <- X.res
        NULL
    })
    x_mat
}


# ---------------------------------------------------------------------------
#' @describeIn anglemania_utils Utility function to extract the genes that
#' have been selected by the anglemania algorithm.
#' @param sce A SingleCellExperiment or SummarizedExperiment object
#' @importFrom S4Vectors metadata
#'
#' @return A character vector of gene names that have been selected by the
#' anglemania algorithm
#' @examples
#' sce <- sce_example()
#' sce <- anglemania(sce, batch_key = "batch")
#' anglemania_genes <- get_anglemania_genes(sce)
#' head(anglemania_genes)
#' length(anglemania_genes)
#' @export
get_anglemania_genes <- function(sce) {
    S4Vectors::metadata(sce)$anglemania$anglemania_genes
}


# ---------------------------------------------------------------------------
#' @describeIn anglemania_utils Utility function to extract the stats of the
#' gene pairs from which the anglemania genes were selected.
#'
#' @param sce A SingleCellExperiment or SummarizedExperiment object
#' @importFrom  S4Vectors metadata
#'
#' @return A data frame of gene pairs from which the anglemania genes
#' were selected
#' @examples
#' sce <- sce_example()
#' sce <- anglemania(sce, batch_key = "batch")
#' anglemania_stats_df <- get_anglemania_stats_df(sce)
#' head(anglemania_stats_df)
#' length(anglemania_stats_df)
#' @export
get_anglemania_stats_df <- function(sce) {
    zscore_mean_threshold <- metadata(sce)$anglemania$params$zscore_mean_threshold
    zscore_sn_threshold <- metadata(sce)$anglemania$params$zscore_sn_threshold
    direction <- metadata(sce)$anglemania$params$direction
    anglemania_genes <- metadata(sce)$anglemania$anglemania_genes
    prefiltered_df <- metadata(sce)$anglemania$prefiltered_df
    prefiltered_df$geneA <-
        S4Vectors::metadata(sce)$anglemania$intersect_genes[
            prefiltered_df$geneA
        ]
    prefiltered_df$geneB <-
        S4Vectors::metadata(sce)$anglemania$intersect_genes[
            prefiltered_df$geneB
        ]
    # Selects the direction of conserved genes
    if (direction == "both") {
        filtered_genes_df <- subset(
            prefiltered_df,
            sn_zscore >= zscore_sn_threshold &
                abs(prefiltered_df$mean_zscore) >= zscore_mean_threshold
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
    # Order data frame
    filtered_genes_df <- filtered_genes_df[
        order(abs(filtered_genes_df$mean_zscore), decreasing = TRUE),
    ]
    filtered_genes_df <- filtered_genes_df |>
        filter(geneA %in% anglemania_genes | geneB %in% anglemania_genes)
    
    return(filtered_genes_df)
}



