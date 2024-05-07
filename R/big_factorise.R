#' Factorise angle matrices
#'
#' @description
#' Generates two a one-hot encoded square matrixces
#' recording the presence/absence (1/0) of sharp and blunt
#' critical angles.
#'
#' @details
#' *factorise* extracts angles between genes, estimates
#' critical angles by approximation of the angle distribution,
#' and records angles passing the critical threshold into
#' a sparse matrix.
#'
#' @param x_mat Matrix. Contains normalised and scaled gene
#'   expression.
#' @param name character. The name of the dataset.
#' @param fdr_threshold double. Fraction of the correlation
#'   to be cut from both sides of an approximated angle
#'   distribution.
#' @return list. First two elements are sparse matrices
#'   containing the significant sharp and blunt angles
#'   between genes. Third and fourth elements are lists
#'   with angles statistics and paths to the values of
#'   critical mangles.
#' @export big_factorise
big_factorise <- function(x_mat, # nolint
                          name) {
    x_mat <- sparse_to_fbm(x_mat)
    # initialize empty FBM with same dimensions as x_mat to store permuted correlation matrix
    x_mat_perm <- bigstatsr::FBM(nrow = nrow(x_mat), ncol = ncol(x_mat))
    
    # permute matrix
    bigstatsr::big_apply(x_mat, a.FUN = function(X, ind) {
        set.seed(1)
        X.sub <- X[, ind, drop = FALSE]
        X.sub <- apply(X.sub, 2, sample)
        x_mat_perm[, ind] <- X.sub
        NULL
    }, a.combine = "c", block.size = 200)

    # compute correlation matrix for both original and permuted matrix
    x_mat_corr <- big_extract_corr(x_mat)
    x_mat_perm_corr <- big_extract_corr(x_mat_perm)

    # get the mean and standard deviation of the permuted correlation matrix 
    # will be used to transform the original correlation matrix into zscores
    dstat <- get_dstat(x_mat_perm_corr)

    # transform original correlation matrix into zscores
    bigstatsr::big_apply(x_mat_corr, a.FUN = function(X, ind) {
        X[, ind] <- (X[, ind, drop = FALSE] - dstat$mean) / dstat$sd
    }, a.combine = "c", block.size = 200)

    return(x_mat_corr)
}


