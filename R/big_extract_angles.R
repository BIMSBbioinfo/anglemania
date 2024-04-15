#'
#' @description
#' Constructs a matrix of angles between gene pairs. Basically this is a converted correlation matrix
#' leverages file-backed matrices from the bigstatsr package.
#'
#' @param x_mat matrix. Contains normalised and scaled gene
#'   expression, where rows are genes and columns are samples.
#' @param gene_names list. A list of vectors specifying the gene_names of the datasets
#' for Seurat object this would be the rownames of the Seurat object (rownames(se))
#' Specified before in \code{\link{big_anglemanise}}.
#' @return FBM (bigstatsr file-backed matrix). Square matrix containing angles between
#' vectors of gene expression (rwos of the input matrix).
#' @export big_extract_angles
big_extract_angles <- function(
    x_mat) {
    # validate inputs
    if (!is.matrix(x_mat)) {
        stop("x_mat has to be a matrix")
    }
    X <- as_FBM(x_mat)
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
    # probably don't need scaling. only do normalization
    # scale_factors <- big_scale(center = TRUE, scale = TRUE)(X)


    # message("Computing cosine distances...")
    # Calculate correlation matrix
    X <- big_transpose(X)

    X <- big_cor(X, block.size = 1000)
    print("inside big_extract_angles")
    print(X[1:5, 1:5])
    # # Convert correlation matrix to an angle matrix
    # NOTE: use correlation for now
    # big_apply(X, a.FUN = function(X, ind) {
    #     # TODO: make something like a progress bar
    #     # access a subset of columns as a standard R matrix
    #     X.sub <- X[, ind, drop = FALSE]

    #     X[, ind] <- round(acos(X.sub) * 180 / pi, 2)
    #     NULL
    # }, a.combine = "c", block.size = 500)
    diag(X) <- NA

    return(X)
}
