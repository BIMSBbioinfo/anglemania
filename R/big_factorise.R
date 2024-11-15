#' Factorize Angle Matrices into Z-Scores
#'
#' @description
#' `big_factorise` computes the angle matrix of the input gene expression
#' matrix using the specified method, performs permutation to create a null
#' distribution, and transforms the correlations into z-scores. This function
#' is optimized for large datasets using the \pkg{bigstatsr} package.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item **Permutation**: The input matrix is permuted column-wise to disrupt
#'     existing angles, creating a null distribution.
#'   \item **Angle Computation**: Computes the angle matrix for both the
#'     original and permuted matrices using \code{\link{big_extract_corr}}.
#'   \item **Method-Specific Processing**:
#'   \itemize{
#'     \item If \code{method = "diem"}, computes Euclidean distances and scales
#'       the angles accordingly, based on the methodology from the DIEM
#'       algorithm (\url{https://bytez.com/docs/arxiv/2407.08623/paper}).
#'     \item For other methods (\code{"pearson"}, \code{"spearman"}),
#'       statistical measures are computed from the permuted data.
#'   }
#'   \item **Statistical Measures**: Calculates mean, variance, and standard
#'     deviation using \code{\link{get_dstat}}.
#'   \item **Z-Score Transformation**: Transforms the original angle matrix into
#'     z-scores.
#' }
#' This process allows for the identification of invariant gene-gene
#' relationships by comparing them to a null distribution derived from the
#' permuted data.
#'
#' @param x_mat A \code{\link[bigstatsr]{FBM}} object representing the
#'   normalized and scaled gene expression matrix.
#' @param method A character string specifying the method for calculating the
#'   relationship between gene pairs. Default is \code{"pearson"}. Other options
#'   include \code{"spearman"} and \code{"diem"} (see
#'   \url{https://bytez.com/docs/arxiv/2407.08623/paper}).
#' @param seed An integer value for setting the seed for reproducibility during
#'   permutation. Default is \code{1}.
#'
#' @return An \code{\link[bigstatsr]{FBM}} object containing the
#'   z-score-transformed angle matrix.
#'
#' @importFrom bigstatsr FBM big_apply
#'
#' @seealso
#' \code{\link{big_extract_corr}},
#' \code{\link{get_dstat}},
#' \code{\link[bigstatsr]{big_apply}},
#' \code{\link[bigstatsr]{FBM}}
#'
#' @export
big_factorise <- function(
    x_mat,
    method = "pearson",
    seed = 1) {
  # Initialize empty FBM to store permuted correlation matrix
  x_mat_perm <- bigstatsr::FBM(
    nrow = nrow(x_mat),
    ncol = ncol(x_mat)
  )

  # Permute matrix
  set.seed(seed)
  bigstatsr::big_apply(
    x_mat,
    a.FUN = function(X, ind) {
      X.sub <- X[, ind, drop = FALSE]
      X.sub <- apply(X.sub, 2, sample)
      x_mat_perm[, ind] <- X.sub
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )

  # Compute correlation matrix for both original and permuted matrix
  x_mat_corr <- big_extract_corr(x_mat, method = method)
  x_mat_perm_corr <- big_extract_corr(x_mat_perm, method = method)

  if (method == "diem") {
    bigstatsr::big_apply(
      x_mat_corr,
      a.FUN = function(X, ind) {
        # 1. Calculate Euclidean distance from Pearson correlation:
        #    dij = sqrt(2 * (1 - rij))
        X.sub <- X[, ind, drop = FALSE]
        X.sub <- sqrt(2 * (1 - X.sub))
        x_mat_corr[, ind] <- X.sub
        NULL
      },
      a.combine = "c",
      block.size = 1000
    )

    dstat <- get_dstat(x_mat_corr)

    scale_factor <- (dstat$min - dstat$max) / dstat$var
    bigstatsr::big_apply(
      x_mat_corr,
      a.FUN = function(X, ind) {
        # 2. Calculate DIEM by subtracting the mean and scaling
        X.sub <- X[, ind, drop = FALSE]
        X.sub <- scale_factor * (X.sub - dstat$mean)
        x_mat_corr[, ind] <- X.sub
        NULL
      },
      a.combine = "c",
      block.size = 1000
    )
  } else {
    dstat <- get_dstat(x_mat_perm_corr)
  }

  # Transform original correlation matrix into z-scores
  bigstatsr::big_apply(
    x_mat_corr,
    a.FUN = function(X, ind) {
      X[, ind] <- (X[, ind, drop = FALSE] - dstat$mean) / dstat$sd
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )

  return(x_mat_corr)
}
