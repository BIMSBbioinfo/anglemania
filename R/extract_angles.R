#' Calculate cosine angle between genes
#' @description
#' Constructs a matrix of gene-gene relationships based on distance
#' metrics.
#' @details
#' The function returns the gene-gene angle matrix as an
#' \code{\link[bigstatsr]{FBM}} object.
#'
#' @param x_mat An \code{\link[bigstatsr]{FBM}} object containing raw gene
#'   expression data, where rows correspond to genes and columns to samples.
#'   The data will be normalized and scaled within the function.
#' @param method A character string specifying the method to compute the
#'   gene-gene relationships. Options are:
#'   \itemize{
#'     \item \code{"cosine"} (default): Computes the cosine angle between
#'       genes.
#'     \item \code{"spearman"}: Computes the Spearman rank correlation
#'       coefficient by rank-transforming the data before computing the
#'       correlation.
#'     \item \code{"diem"}: Computes the Dimension Insensitive Euclidean
#'       Metric between genes. Note that this is done in the
#'       \code{\link{factorise}} function.
#'   }
#'
#' @return An \code{\link[bigstatsr]{FBM}} object containing the gene-gene
#'   correlation matrix. The matrix is square with dimensions equal to the
#'   number of genes and contains the pairwise correlations between genes.
#'   The diagonal elements are set to \code{NA}.
#'
#' @importFrom bigstatsr FBM big_apply big_transpose big_cor
#' @examples 
#'
#' mat <- matrix(
#'  c(
#'      5, 3, 0, 0,
#'      0, 0, 0, 3,
#'      2, 1, 3, 4,
#'      0, 0, 1, 0,
#'      1, 2, 1, 2,
#'      3, 4, 3, 4
#'    ),
#'    nrow = 6, # 6 genes
#'    ncol = 4, # 4 cells
#'    byrow = TRUE
#' )
#'
#' mat <- bigstatsr::FBM(nrow = nrow(mat), ncol = ncol(mat), init = mat)
#'
#' angle_mat <- extract_angles(mat)
#' angle_mat[]
#' @seealso
#' \code{\link[bigstatsr]{big_apply}},
#' \code{\link[bigstatsr]{big_cor}},
#' \code{\link[bigstatsr]{FBM}},
#' \code{\link{factorise}}
#'
#' @export
extract_angles <- function(
    x_mat,
    method = "cosine") {
  bigstatsr::big_apply(x_mat, a.FUN = function(X, ind) {
    X.sub <- X[, ind, drop = FALSE]
    # Normalize the data:
    #   divide gene counts by the number of total counts per cell,
    #   multiply by 10,000 (scaling factor like in Seurat)
    X.sub <- t(t(X.sub) / colSums(X.sub) * 10000)
    X.sub <- log1p(X.sub)

    X[, ind] <- X.sub
    NULL
  })

  # First transpose the matrix because big_cor calculates the covariance
  # (X^T X)
  x_mat <- bigstatsr::big_transpose(x_mat)

  # Transform to ranks if method is spearman
  if (method == "spearman") {
    bigstatsr::big_apply(x_mat, a.FUN = function(X, ind) {
      X.sub <- X[, ind, drop = FALSE]
      X.sub <- apply(X.sub, 2, rank)
      x_mat[, ind] <- X.sub
      NULL
    })
  }

  x_mat <- bigstatsr::big_cor(x_mat, block.size = 1000)
  # The big_cor function from bigstatsr scales and centers the
  # count matrix and calculates the covariance (cross product X^T X)
  diag(x_mat) <- NA

  return(x_mat)
}
