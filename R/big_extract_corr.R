#'
#' @description
#' Constructs a matrix of angles between gene pairs. Basically this is a converted correlation matrix
#' leverages file-backed matrices from the bigstatsr package.
#'
#' @param x_mat matrix. Contains normalised and scaled gene
#'   expression, where rows are genes and columns are samples.
#' @param gene_names vector. A character vector specifying the gene_names of the dataset
#' for Seurat object this would be the rownames of the Seurat object (rownames(se)) 
#' ==> when using big_factorise, those are already extracted
#' Specified before in \code{\link{big_factorise}}.
#' @return FBM (bigstatsr file-backed matrix). Square matrix containing angles between
#' vectors of gene expression (rwos of the input matrix).
#' @export big_extract_angles
big_extract_corr <- function(
    x_mat) {
    # validate inputs
    if (!is.matrix(x_mat)) {
        stop("x_mat has to be a matrix")
    }
    X <- as_FBM(x_mat)
    df <- ncol(X)-2 # important for the degrees of freedom when computing the p-values
    print(paste0("computing p-values for ", n, " cells"))
    # log normalize the data
    big_apply(X, a.FUN = function(X, ind) {
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
    X <- big_transpose(X)

    X <- big_cor(X, block.size = 1000)
    # the big_cor function from bigstatsr basically scales and centers the count matrix and calculates the covariance (cross product XT*X)
    diag(X) <- NA

    # Calculate p-values and adjust them using the Benjamini-Hochberg method
    big_apply(X, a.FUN = function(X, ind) {
        X.sub <- X[, ind, drop = FALSE]
        computePValues(X.sub, df = df)
        # a C++ function that computes the p-values from the correlation matrix and replaces the values of the matrix in place.
        # uses the CDF of a normal distribution not a students-t distribution. Justified by the fact that with higher df towards 30, students t CDF is very close to normal distribution.
        X.sub <- p.adjust(X.sub, method = "BH")
        X[, ind] <- X.sub
        NULL
    }, a.combine = "c", block.size = 1000)

    return(X)
}
