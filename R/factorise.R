#' Factorize Angle Matrices into Z-Scores
#'
#' @description
#' `factorise` computes the angle matrix of the input gene expression
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
#'     original and permuted matrices using \code{\link{extract_angles}}.
#'   \item **Method-Specific Processing**:
#'   \itemize{
#'     \item For other methods (\code{"cosine"}, \code{"spearman"}),
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
#'   relationship between gene pairs. Default is \code{"cosine"}. Other options
#'   include \code{"spearman"}
#' @param seed An integer value for setting the seed for reproducibility during
#'   permutation. Default is \code{1}.
#' @param permute_row_or_column Character "row" or "column", whether permutations should be executed row-wise or column wise. Default is \code{"column"}
#' @param permutation_function Character "sample" or "permute_nonzero". If sample, then sample is used for constructing background distributions. If permute_nonzero, then only non-zero values are permuted. Default is \code{"sample"}
#' @param normalization_method Character "divide_by_total_counts" or
#'   "scale_by_total_counts". Default is \code{"divide_by_total_counts"}
#' @return An \code{\link[bigstatsr]{FBM}} object containing the
#'   z-score-transformed angle matrix.
#'
#' @importFrom bigstatsr FBM big_apply rows_along
#' @importFrom checkmate assertClass assertString assertChoice
#' @importFrom withr with_seed
#' @importFrom digest digest
#'
#' @examples
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
#' # Run factorise with method "cosine" and a fixed seed
#' result_fbm <- factorise(mat, method = "cosine", seed = 1)
#' result_fbm[]
#' @seealso
#' \code{\link{extract_angles}},
#' \code{\link{get_dstat}},
#' \code{\link[bigstatsr]{big_apply}},
#' \code{\link[bigstatsr]{FBM}}
#'
#' @export
factorise <- function(
    x_mat,
    method = "cosine",
    seed = 1,
    permute_row_or_column = "column",
    permutation_function = "sample",
    normalization_method = "divide_by_total_counts"
  ) {
  # Validate input
  checkmate::assertClass(x_mat, "FBM")
  checkmate::assertString(method)
  checkmate::assertChoice(method, c("cosine", "spearman"))
  checkmate::assertString(permute_row_or_column)
  checkmate::assertChoice(permute_row_or_column, c("row", "column"))
  checkmate::assertString(permutation_function)
  checkmate::assertChoice(permutation_function, c("sample", "permute_nonzero"))
  tmpfile = tempfile()
  x_mat_perm <- bigstatsr::FBM(
    nrow = nrow(x_mat),
    ncol = ncol(x_mat),
    backingfile = file.path(tmpfile, digest::digest(tmpfile, length = 10))
  )

  if(permutation_function == "sample"){
    permutation_function = base::sample
  }else{
    permutation_function = permute_nonzero
  }

  # normalizes the matrix before permutation
  x_mat = normalize_matrix(x_mat, normalization_method = normalization_method)

  # default permutation is by columns
  ind_fun = bigstatsr::cols_along(x_mat)
  a_fun = function(X, ind) {  
      X.sub <- X[, ind, drop = FALSE]
      X.sub <- apply(X.sub, 2, permutation_function)
      x_mat_perm[, ind] <- X.sub
      NULL
    }

  # permute by rows
  if(permute_row_or_column == "row"){
    ind_fun = bigstatsr::rows_along(x_mat)
    a_fun = function(X, ind){
      X.sub <- X[ind, ,drop = FALSE]
      X.sub <- t(apply(X.sub, 1, permutation_function))
      x_mat_perm[ind, ] <- X.sub
    }
  }
  withr::with_seed(seed, {
    # Permute matrix
    bigstatsr::big_apply(
      x_mat,
      a.FUN = a_fun,
      a.combine = "c",
      ind = ind_fun,
      block.size = 1000
    )
  })

  # Compute correlation matrix for both original and permuted matrix
  x_mat_corr <- extract_angles(x_mat, method = method)
  x_mat_perm_corr <- extract_angles(x_mat_perm, method = method)

  # Transform original correlation matrix into z-scores.
  dstat <- get_dstat(x_mat_perm_corr)
  bigstatsr::big_apply(
    x_mat_corr,
    a.FUN = function(X, ind) {
      zscores <- (X[, ind, drop = FALSE] - dstat$mean[ind]) / dstat$sd[ind]
      # this is needed because sometimes all the correlation matrices are 0, and then
      # the mean is zero
      zscores[is.na(zscores)] = 0
      X[, ind] <- zscores
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )
  return(x_mat_corr)
}


# ---------------------------------------------------------------------------
#' permute non-zero elements of vectors
#' @description
#' Permutes the non-zero elements of a numeric vector.
#' @param v A numeric vector.
#' @return A numeric vector with non-zero elements permuted.
#' @examples
#' v <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) # put in some zeroes
#' v = c(v, 0, 0, 0)
#' permute_nonzero(v)
#' @export
permute_nonzero = function(v) {
  ind = v > 0
  v[ind] = sample(v[ind])
  v
}