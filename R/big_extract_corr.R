#'
#' @description
#' Constructs a matrix of correlations between gene pairs.
#' 
#' @param x_mat matrix. Contains normalised and scaled gene
#'   expression, where rows are genes and columns are samples.
#' @param method character. Method to be used for calculating the relationship of a gene pair.
#'   Default is pearson correlation. Other options are 'pearson' and 'spearman' (or 'diem') but 
#'      the code only changes for spearman.
#' @return FBM (bigstatsr file-backed matrix). Square matrix containing correlations between
#' vectors of gene expression (rwos of the input matrix).
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

    # Calculate correlation matrix
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
