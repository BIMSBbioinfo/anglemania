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
    method = "cosine"
  ) {
  checkmate::assertChoice(method, c("cosine", "spearman"))
  checkmate::assertClass(x_mat, "FBM")
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

  # x_mat <- bigstatsr::big_cor(x_mat, block.size = 1000)
  x_mat <- big_cor_no_warning(x_mat, block.size = 1000)
  # x_mat <- angl_cor(x_mat, block.size = 1000)
  # The big_cor function from bigstatsr scales and centers the
  # count matrix and calculates the covariance (cross product X^T X)
  diag(x_mat) <- NA

  # replaces NaN with NA values in the matrix
  bigstatsr::big_apply(x_mat, a.FUN = function(X, ind) {
    X.sub <- apply(X[, ind, drop = FALSE], 2, function(x) {
      xx = x
      xx[is.nan(xx)] = NA
      xx
    })
    x_mat[, ind] <- X.sub
    NULL
  })

  return(x_mat)
}
