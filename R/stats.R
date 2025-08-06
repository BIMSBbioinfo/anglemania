# ---------------------------------------------------------------------------
# Statistics functions for the anglemania package
# ---------------------------------------------------------------------------
#' @title Statistics functions for the anglemania package
#' @description A collection of utility functions used within the
#' \pkg{anglemania} package for calculating statistics based
#' on the results created during the anglemania function.
#' @name statistics
#' @rdname statistics
#' @keywords internal
NULL
# ---------------------------------------------------------------------------
#' @describeIn statistics Compute mean, standard deviation, variance,
#' min, and max of a correlation matrix stored as an
#' \code{\link[bigstatsr]{FBM}}.
#' @param corr_matrix An \code{\link[bigstatsr]{FBM}} object.
#' @return A list with statistical measures including \code{mean}, \code{sd},
#'   \code{var}, \code{min}, and \code{max}.
#' @importFrom bigstatsr big_apply
#' @examples
#' s_mat <- Matrix::rsparsematrix(nrow = 10, ncol = 5, density = 0.3)
#' fbm_mat <- sparse_to_fbm(s_mat)
#' result <- get_dstat(fbm_mat)
#' str(result)
#' result
#' @seealso \code{\link[bigstatsr]{big_apply}}, \code{\link[bigstatsr]{FBM}}
#' @export
get_dstat <- function(corr_matrix) {
  # Check if FBM is used as input
  if (!inherits(corr_matrix, "FBM")) {
    stop("corr_matrix has to be a FBM object")
  }
  stats_funs <- c("mean", "var", "min", "max")
  # min/max give warnings if value NA -> suppressWarnings
  dstat <- suppressWarnings({
    lapply(setNames(stats_funs, stats_funs), function(fun) {
      stat_value <- bigstatsr::big_apply(
        corr_matrix,
        a.FUN = function(X, ind) {
          apply(X[, ind, drop = FALSE], 2, fun, na.rm = TRUE)
        },
        a.combine = "c",
        block.size = 1000
      )
      stat_value <- replace_with_na(stat_value)
      stat_value
    })
  })

  # Calculate standard deviation
  dstat$sd <- sqrt(dstat$var)

  return(dstat)
}

# ---------------------------------------------------------------------------
#' @describeIn statistics Calculates the element-wise mean from a list
#' of \code{bigstatsr::FBM} objects.
#'
#' @param matrix_list A list of \code{bigstatsr::FBM} objects.
#'   In this case, the FBMs are the angle matrices computed in \code{factorise}.
#' @param weights A numeric vector of weights for each dataset or batch.
#' @return A new \code{bigstatsr::FBM} object containing the mean values.
#' @importFrom bigstatsr FBM
#' @examples
#' # Create FBMs
#' mat1 <- matrix(1:9, nrow = 3)
#' mat2 <- matrix(1:9, nrow = 3)
#'
#' fbm1 <- bigstatsr::FBM(nrow = nrow(mat1), ncol = ncol(mat1), init = mat1)
#' fbm2 <- bigstatsr::FBM(nrow = nrow(mat2), ncol = ncol(mat2), init = mat2)
#'
#' # Create weights
#' weights <- c(batch1 = 0.5, batch2 = 0.5)
#'
#' # Create the list of FBMs
#' fbm_list <- list(batch1 = fbm1, batch2 = fbm2)
#'
#' big_mat_list_mean(fbm_list, weights)
#' @export
big_mat_list_mean <- function(matrix_list, weights, verbose = TRUE) {
  # Check if all matrices in the list have the same dimensions
  if (
    !all(sapply(
      matrix_list,
      function(x) {
        identical(dim(matrix_list[[1]]), dim(x))
      }
    ))
  ) {
    stop("All matrices in the list must have the same dimensions.")
  }

  n_col <- ncol(matrix_list[[1]])
  n_row <- nrow(matrix_list[[1]])
  mat_mean_zscore <- bigstatsr::FBM(n_row, n_col)

  bigstatsr::big_apply(
    mat_mean_zscore,
    a.FUN = function(X, ind) {
      X.sub <- X[, ind, drop = FALSE]
      wrap_mean <- function(final_mat, batch) {
        batch_mat <- matrix_list[[batch]][,
          ind,
          drop = FALSE
        ]
        # calculates the number of samples in which the feature is
        # present, weighted by the dataset weight
        nmat <- (!is.na(batch_mat) + 0) * weights[batch]
        final_mat <- final_mat + batch_mat
        final_mat[is.na(final_mat)] <- 0
        list(final_mat, nmat)
      }
      # Run function in Reduce statement on names of weights vector
      lmats <- lapply(names(weights), function(batch) {
        wrap_mean(X.sub, batch)
      })
      # this sums up the z-scores across samples
      m_sum <- Reduce("+", lapply(lmats, '[[', 1))
      # this gets the number of samples in which the feature was present
      m_n <- Reduce("+", lapply(lmats, '[[', 2))
      # Divide by zero shouldn't happen because the features have already been filtered properly
      X[, ind] <- m_sum / m_n # Already weighted, no need to divide
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )

  return(mat_mean_zscore)
}

# ---------------------------------------------------------------------------
#' @describeIn statistics Calculate mean, standard deviation, and SNR
#' across a list of FBMs stored in an \code{anglemania_object}.
#'
#' @param matrix_list A list of \code{bigstatsr::FBM} objects.
#' @param weights A numeric vector of weights for each dataset or batch.
#' @return A list containing three matrices: \code{mean_zscore},
#'   \code{sds_zscore}, and \code{sn_zscore} which are later used to
#' filter gene pairs based on the absolute mean z-score and signal-to-noise ratio
#' of the angles.
#' @importFrom bigstatsr FBM big_apply
#' @examples
#' library(SingleCellExperiment)
#' library(S4Vectors)
#' sce <- sce_example()
#' sce <- anglemania(sce, batch_key = "batch")
#' matrix_list <- metadata(sce)$anglemania$matrix_list
#' weights <- setNames(
#'   S4Vectors::metadata(sce)$anglemania$weights$weight,
#'   S4Vectors::metadata(sce)$anglemania$weights$anglemania_batch
#' )
#' list_stats <- get_list_stats(matrix_list, weights)
#' names(list_stats)
#' list_stats$mean_zscore[1:5, 1:5]
#' @seealso \code{\link[bigstatsr]{big_apply}}, \code{\link[bigstatsr]{FBM}}
#' @export
get_list_stats <- function(matrix_list, weights, verbose = TRUE) {
  # Check if all matrices in the list have the same dimensions
  if (
    !all(sapply(
      matrix_list,
      function(x) {
        identical(dim(matrix_list[[1]]), dim(x))
      }
    ))
  ) {
    stop("All matrices in the list need to have the same dimensions.")
  }

  vmessage(verbose, "Weighting matrix_list...")
  vmessage(verbose, "Calculating mean...")
  mat_mean_zscore <- big_mat_list_mean(matrix_list, weights)

  n_col <- ncol(matrix_list[[1]])
  n_row <- nrow(matrix_list[[1]])
  mat_sds_zscore <- bigstatsr::FBM(n_row, n_col)

  vmessage(verbose, "Calculating sds...")
  # calculates the weighted sds for the big matrix
  bigstatsr::big_apply(
    mat_sds_zscore,
    a.FUN = function(X, ind) {
      X.sub <- X[, ind, drop = FALSE]
      wrap_mean <- function(final_mat, batch) {
        batch_mat <- matrix_list[[batch]][,
          ind,
          drop = FALSE
        ]

        # calculates the number of samples in which the feature is present, weighted
        # by the dataset weight
        nmat <- (!is.na(batch_mat) + 0) * weights[batch]

        # calculates the weighted standard deviation
        final_mat <- final_mat +
          (batch_mat - mat_mean_zscore[, ind, drop = FALSE])^2 *
            weights[batch]
        final_mat[is.na(final_mat)] <- 0
        list(final_mat, nmat)
      }
      # Run function in Reduce statement on names of weights vector
      lmats <- lapply(names(weights), function(batch) {
        wrap_mean(X.sub, batch)
      })
      # this sums up the z-scores across samples
      m_sd <- Reduce("+", lapply(lmats, '[[', 1))
      # this gets the number of samples in which the feature was present
      w_m_sum <- Reduce("+", lapply(lmats, function(x) x[[2]]))
      w_m_sum_sq <- Reduce("+", lapply(lmats, function(x) x[[2]]^2))

      # creates an unbiased esitmator for a sample based wighted variance calculation
      denominator <- w_m_sum - (w_m_sum_sq / w_m_sum)

      # Divide by zero shouldn't happen because the features have already been filtered properly
      X[, ind] <- sqrt(m_sd / denominator) # Already weighted, no need to divide
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )

  mat_sn_zscore <- bigstatsr::FBM(n_row, n_col)
  diag(mat_sds_zscore) <- NA # Avoid dividing by zero
  bigstatsr::big_apply(
    mat_sn_zscore,
    a.FUN = function(X, ind) {
      X[, ind] <- abs(
        mat_mean_zscore[, ind, drop = FALSE]
      ) /
        mat_sds_zscore[, ind, drop = FALSE]
      X[, ind] <- replace_with_na(X[, ind])
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )

  res <- list(
    mean_zscore = mat_mean_zscore,
    sds_zscore = mat_sds_zscore,
    sn_zscore = mat_sn_zscore
  )

  return(res)
}
