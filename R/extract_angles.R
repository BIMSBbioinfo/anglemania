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
#' @param n_threads integer. Number of cores to use for
#'   the computation of the cosine distances.
#' @return matrix. Square matrix containing angles between
#' vectors of gene expression (rwos of the input matrix).
#' @export extract_angles
extract_angles <- function(x_mat, #nolint
                           n_threads = 16) {
  x_mat <- as.matrix(x_mat)
  x_mat_ang <- round(
                     acos(
                       1 -
                         parallelDist::parDist(x_mat,
                                               method = "cosine",
                                               diag = FALSE,
                                               upper = FALSE,
                                               threads = n_threads)
                     ) / pi * 180,
                     2)
  x_mat_ang <- dist2mat(x_mat_ang, 256) # efficient conversion
  return(x_mat_ang)
}