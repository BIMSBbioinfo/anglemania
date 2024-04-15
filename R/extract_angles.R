#' Construct a matrix of angles between genes
#'
#' @description
#' Constructs a matrix of angles between genes from an input
#' gene expression matrix.
#'
#' @details
#' This is a **low-level** function that utilises functionality
#' of **parallelDist** pacakge to compute cosine distances
#' between vectors of gene expression. Since **parallelDist**´s#
#' outputs a dist object, an Rcpp functions is used to convert it
#' back to matrix.
#' @importFrom parallelDist parDist
#' @importFrom Rcpp evalCpp
#' @param x_mat matrix. Contains normalised and scaled gene
#'   expression, where rows are genes and columns are samples.
#' param n_cores integer. Number of cores to use for
#'   the computation of the cosine distances.
#' @return matrix. Square matrix containing angles between
#' vectors of gene expression (rwos of the input matrix).
#' @export extract_angles
extract_angles <- function(x_mat # nolint
) {
  x_mat <- as.matrix(x_mat)
  x_mat_ang <- round(
    acos( # acos returns values between 0 and pi ==> convert to degrees later
      1 - # 1 - cosine distance = cosine similarity
        parallelDist::parDist(x_mat, # compute cosine distance
          method = "cosine",
          diag = FALSE,
          upper = FALSE,
          threads = 1
        )
    ) / pi * 180, # convert to degrees
    2
  )
  x_mat_ang <- dist2mat(x_mat_ang, 256) # efficient conversion
  return(x_mat_ang)
}
