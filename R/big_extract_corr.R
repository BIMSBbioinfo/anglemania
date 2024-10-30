#'
#' @description
#' Constructs a matrix of gene-gene relationships based on distance metrics.
#' @details
#' The function returns the gene-gene angle matrix as an \code{\link[bigstatsr]{FBM}} object.
#'
#' @param x_mat An \code{\link[bigstatsr]{FBM}} object containing raw gene expression data, where rows correspond to genes and columns to samples. The data will be normalized and scaled within the function.
#' @param method A character string specifying the method to compute the  Options are:
#'   \itemize{
#'     \item \code{"pearson"} (default): Computes the cosine angle between genes.
#'     \item \code{"spearman"}: Computes the Spearman rank correlation coefficient by rank-transforming the data before computing the correlation.
#'     \item \code{"diem"}: computes the Dimension Insensitive Euclidean Metric between genes. Note that this is done in the \code{\link{big_factorise}} function.
#'   }
#'
#' @return An \code{\link[bigstatsr]{FBM}} object containing the gene-gene correlation matrix. The matrix is square with dimensions equal to the number of genes and contains the pairwise correlations between genes. The diagonal elements are set to \code{NA}.
#'
#' @importFrom bigstatsr FBM big_apply big_transpose big_cor
#'
#' @seealso
#' \code{\link[bigstatsr]{big_apply}}, \code{\link[bigstatsr]{big_cor}}, \code{\link[bigstatsr]{FBM}}, \code{\link{big_factorise}}
#'
#' @examples
#' \dontrun{
#' # Load necessary library
#' library(bigstatsr)
#'
#' # Create a random gene expression FBM for demonstration
#' n_genes <- 1000
#' n_samples <- 500
#' set.seed(123)
#' x_mat <- FBM(n_genes, n_samples, init = rpois(n_genes * n_samples, lambda = 5))
#'
#' # Compute the gene-gene correlation matrix using Pearson correlation
#' corr_matrix <- big_extract_corr(x_mat, method = "pearson")
#'
#' # Compute the gene-gene correlation matrix using Spearman correlation
#' corr_matrix_spearman <- big_extract_corr(x_mat, method = "spearman")
#'
#' # Access a subset of the correlation matrix
#' corr_subset <- corr_matrix[1:5, 1:5]
#' print(corr_subset)
#' }
#'
#' @export big_extract_corr
big_extract_corr <- function(
    x_mat,
    method = "pearson") {
        
    bigstatsr::big_apply(x_mat, a.FUN = function(X, ind) {
        X.sub <- X[, ind, drop = FALSE]
        # normalize the data:
        #   divide gene counts by the number of total counts per cell,
        #   multiply by 10,000 (scaling factor like in Seurat)
        X.sub <- t(t(X.sub) / colSums(X.sub) * 10000)
        X.sub <- log1p(X.sub)

        X[, ind] <- X.sub
        NULL
    })

    # first transpose the matrix because big_cor calculates the covariance (XT*X)
    x_mat <- bigstatsr::big_transpose(x_mat)

    # transform to ranks if method is spearman
    if (method == "spearman"){
        big_apply(x_mat, a.FUN = function(X, ind) {
            X.sub <- X[, ind, drop = FALSE]
            X.sub <- apply(X.sub, 2, rank)
            x_mat[, ind] <- X.sub
            NULL
        })
    }
    x_mat <- bigstatsr::big_cor(x_mat, block.size = 1000)
    # the big_cor function from bigstatsr basically scales and centers the count matrix and calculates the covariance (cross product XT*X)
    diag(x_mat) <- NA

    return(x_mat)
}
