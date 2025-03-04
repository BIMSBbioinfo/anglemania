# ---------------------------------------------------------------------------
# Utility functions for the anglemania package
# ---------------------------------------------------------------------------

#' Convert a sparse matrix into a file-backed matrix
#' (\code{\link[bigstatsr]{FBM}})
#'
#' Converts a sparse matrix into an \code{\link[bigstatsr]{FBM}} with efficient
#' memory usage.
#'
#' @param s_mat A sparse matrix.
#' @return An \code{\link[bigstatsr]{FBM}} object from the \pkg{bigstatsr}
#'  package.
#' @importFrom bigstatsr FBM
#' @examples
#' s_mat <- Matrix::rsparsematrix(nrow = 10, ncol = 5, density = 0.3)
#' # Convert the sparse matrix to an FBM using your function
#' fbm_mat <- sparse_to_fbm(s_mat)
#' fbm_mat
#' @export
sparse_to_fbm <- function(s_mat) {
  n <- nrow(s_mat)
  p <- ncol(s_mat)
  X <- bigstatsr::FBM(n, p)

  bigstatsr::big_apply(
    X,
    a.FUN = function(X, ind) {
      X[, ind] <- s_mat[, ind] %>% as.matrix()
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )

  return(X)
}

# ---------------------------------------------------------------------------
#' Compute mean and standard deviation of the correlation matrix
#'
#' Computes the mean and standard deviation of the correlation matrix using the
#' \code{big_apply} function.
#'
#' @param corr_matrix An \code{\link[bigstatsr]{FBM}} object.
#' @return A list with statistical measures including \code{mean}, \code{sd},
#'   \code{var}, \code{sn}, \code{min}, and \code{max}.
#' @importFrom bigstatsr big_apply
#' @examples
#' s_mat <- Matrix::rsparsematrix(nrow = 10, ncol = 5, density = 0.3)
#' # Convert the sparse matrix to an FBM using your function
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

  # Calculate number of non-NA entries
  n <- bigstatsr::big_apply(
    corr_matrix,
    a.FUN = function(X, ind) {
      sum(!is.na(X[, ind, drop = FALSE]))
    },
    a.combine = "sum",
    block.size = 1000
  )

  # Calculate mean
  mean_value <- bigstatsr::big_apply(
    corr_matrix,
    a.FUN = function(X, ind) {
      sum(X[, ind, drop = FALSE], na.rm = TRUE)
    },
    a.combine = "sum",
    block.size = 1000
  ) / n

  # Calculate variance
  var_value <- bigstatsr::big_apply(
    corr_matrix,
    a.FUN = function(X, ind) {
      sum(
        (X[, ind, drop = FALSE] - mean_value)^2,
        na.rm = TRUE
      ) / (n - 2) # minus 2 because it's a symmetric matrix!!
    },
    a.combine = "sum",
    block.size = 1000
  )

  # Calculate standard deviation
  sd_value <- sqrt(var_value)

  # Calculate minimum
  min_value <- bigstatsr::big_apply(
    corr_matrix,
    a.FUN = function(X, ind) {
      min(X[, ind, drop = FALSE], na.rm = TRUE)
    },
    a.combine = "min",
    block.size = 1000
  )

  # Calculate maximum
  max_value <- bigstatsr::big_apply(
    corr_matrix,
    a.FUN = function(X, ind) {
      max(X[, ind, drop = FALSE], na.rm = TRUE)
    },
    a.combine = "max",
    block.size = 1000
  )

  dstat <- list(
    mean = mean_value,
    var  = var_value,
    sd   = sd_value,
    sn   = sd_value / mean_value, # signal-to-noise ratio
    min  = min_value,
    max  = max_value
  )

  return(dstat)
}

# ---------------------------------------------------------------------------
#' Calculates the mean for each element of a matrix from a list of FBMs
#'
#' This function takes an \code{anglemania_object} containing a list of
#' \code{\link[bigstatsr]{FBM}}s and calculates the mean of every element.
#' If the list is empty or the FBMs have different dimensions,
#' it throws an error.
#'
#' @param angl An \code{anglemania_object} containing the list of FBMs.
#'  In this case, the FBMs are the angle matrices computed in \code{factorise}.
#' @return A new \code{\link[bigstatsr]{FBM}} object containing the mean values.
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
#' # Construct the anglemania_object
#' angl <- new(
#'   "anglemania_object",
#'   weights = weights,
#'   matrix_list = fbm_list
#' )
#' big_mat_list_mean(angl)
#' @export
big_mat_list_mean <- function(angl) {
  if (!inherits(angl, "anglemania_object")) {
    stop("angl needs to be an anglemania_object")
  }

  # Check if all matrices in the list have the same dimensions
  if (!all(sapply(
    matrix_list(angl),
    function(x) {
      identical(dim(matrix_list(angl)[[1]]), dim(x))
    }
  ))) {
    stop("All matrices in the list must have the same dimensions.")
  }

  n_col <- ncol(matrix_list(angl)[[1]])
  n_row <- nrow(matrix_list(angl)[[1]])
  mat_mean_zscore <- bigstatsr::FBM(n_row, n_col)

  bigstatsr::big_apply(
    mat_mean_zscore,
    a.FUN = function(X, ind) {
      X.sub <- X[, ind, drop = FALSE]
      wrap_mean <- function(final_mat, batch) {
        batch_mat <- matrix_list(angl)[[batch]][
          , ind,
          drop = FALSE
        ] * angl@weights[batch]
        final_mat <- final_mat + batch_mat
      }
      # Run function in Reduce statement on names of weights vector
      X.sub <- Reduce(
        wrap_mean,
        names(angl@weights),
        init = X.sub
      )
      X[, ind] <- X.sub # Already weighted, no need to divide
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )

  return(mat_mean_zscore)
}

# ---------------------------------------------------------------------------
#' Calculate statistical measures from a list of FBMs
#'
#' Computes the mean, standard deviations, and signal-to-noise ratio (SNR)
#' for each element across a list of FBMs in an \code{anglemania_object}.
#'
#' @param angl An \code{anglemania_object} containing the list of FBMs.
#' @return A list containing three matrices: \code{mean_zscore},
#'   \code{sds_zscore}, and \code{sn_zscore}.
#' @importFrom bigstatsr FBM big_apply
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl <- anglemania(angl)
#' list_stats(angl) <- get_list_stats(angl)
#' str(list_stats(angl))
#' @seealso \code{\link[bigstatsr]{big_apply}}, \code{\link[bigstatsr]{FBM}}
#' @export
get_list_stats <- function(angl) {
  if (!inherits(angl, "anglemania_object")) {
    stop("angl needs to be an anglemania_object")
  }

  # Check if all matrices in the list have the same dimensions
  if (!all(sapply(
    matrix_list(angl),
    function(x) {
      identical(dim(matrix_list(angl)[[1]]), dim(x))
    }
  ))) {
    stop("All matrices in the list need to have the same dimensions.")
  }

  message("Weighting matrix_list...")
  message("Calculating mean...")

  mat_mean_zscore <- big_mat_list_mean(angl)

  n_col <- ncol(matrix_list(angl)[[1]])
  n_row <- nrow(matrix_list(angl)[[1]])
  mat_sds_zscore <- bigstatsr::FBM(n_row, n_col)

  message("Calculating sds...")

  bigstatsr::big_apply(
    mat_sds_zscore,
    a.FUN = function(X, ind) {
      wrap_sds <- function(final_mat, batch) {
        batch_mat <- matrix_list(angl)[[batch]]
        final_mat <- final_mat + (
          batch_mat[, ind, drop = FALSE] -
            mat_mean_zscore[, ind, drop = FALSE]
        )^2 * angl@weights[batch]
      }
      X.sub <- X[, ind, drop = FALSE]
      X.sub <- Reduce(
        wrap_sds,
        names(angl@weights),
        init = X.sub
      )
      X[, ind] <- sqrt(X.sub)
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
      ) / mat_sds_zscore[, ind, drop = FALSE]
      NULL
    },
    a.combine = "c",
    block.size = 1000
  )

  res <- list(
    mean_zscore = mat_mean_zscore,
    sds_zscore  = mat_sds_zscore,
    sn_zscore   = mat_sn_zscore
  )

  return(res)
}

# ---------------------------------------------------------------------------
#' Extract Unique Genes from Gene Pairs Data Frame
#'
#' Extracts unique gene identifiers from a data frame containing gene pairs
#' and returns a vector of genes up to a specified maximum number.
#'
#' @param dt A data frame containing gene pairs, with columns \code{geneA}
#'   and \code{geneB}.
#' @param max_n_genes An integer specifying the maximum number of unique genes
#'   to return.
#' @return A vector of unique gene identifiers.
#' @details
#' The function combines the \code{geneA} and \code{geneB} columns, extracts
#' unique gene names, and returns the first \code{max_n_genes} genes. If
#' \code{max_n_genes} exceeds the number of unique genes available, all unique
#' genes are returned.
#' @examples
#' gene_pairs <- data.frame(
#'   geneA = c("Gene1", "Gene2", "Gene3", "Gene4"),
#'   geneB = c("Gene3", "Gene4", "Gene5", "Gene6")
#' )
#' unique_genes <- extract_rows_for_unique_genes(
#'   gene_pairs,
#'   max_n_genes = 3
#' )
#' print(unique_genes)
#' @seealso \code{\link{select_genes}}
#' @export
extract_rows_for_unique_genes <- function(dt, max_n_genes) {
  unique_genes <- unique(as.vector(rbind(dt$geneA, dt$geneB)))
  max_genes <- ifelse(
    max_n_genes > length(unique_genes),
    length(unique_genes),
    max_n_genes
  )
  unique_genes <- unique_genes[1:max_genes]
  return(unique_genes)
}


# ---------------------------------------------------------------------------
#' Preselect Genes from an anglemania_object, to make subsequent filtering easier
#'
#' @param angl An \code{anglemania_object} containing statistical
#'   matrices such as mean z-scores and SNR z-scores.
#' @param zscore_mean_threshold Numeric value specifying the threshold for the
#'   absolute mean z-score. Default is 1.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   SNR z-score. Default is 1.
#' @return The input \code{anglemania_object} with the
#'   \code{integration_genes} slot updated to include the selected genes and
#'   their statistical information.
#' @details
#' The function performs the following steps:
#' \enumerate{
#'  \item Checks if the input object is of class \code{\link{anglemania_object-class}}.
#'  \item Identifies gene pairs where both the mean z-score and SNR z-score
#'     exceed the specified thresholds.
#' }
#' @useDynLib anglemania, .registration = TRUE
#' @keywords internal
prefilter_angl <- function(
    angl,
    zscore_mean_threshold = 1,
    zscore_sn_threshold = 1) {
  if (!inherits(angl, "anglemania_object")) {
    stop("angl needs to be an anglemania_object")
  }
  if (zscore_mean_threshold <= 0 || zscore_sn_threshold <= 0) {
    stop("zscore_mean_threshold and zscore_sn_threshold need to be positive")
  }
  prefiltered_dt <- select_genes_cpp(
    BM_sn = list_stats(angl)$sn_zscore,
    BM_mean = list_stats(angl)$mean_zscore,
    zscore_mean_threshold = zscore_mean_threshold,
    zscore_sn_threshold = zscore_sn_threshold
  )
  while (nrow(prefiltered_dt) == 0) {
    if (zscore_mean_threshold <= 0 || zscore_sn_threshold <= 0) {
      stop("zscore_mean_threshold and zscore_sn_threshold need to be positive")
    }
    message("No genes passed the cutoff. Decreasing thresholds by 0.1...")
    zscore_mean_threshold <- zscore_mean_threshold - 0.1
    zscore_sn_threshold <- zscore_sn_threshold - 0.1
    prefiltered_dt <- select_genes_cpp(
      BM_sn = list_stats(angl)$sn_zscore,
      BM_mean = list_stats(angl)$mean_zscore,
      zscore_mean_threshold = zscore_mean_threshold,
      zscore_sn_threshold = zscore_sn_threshold
    )
  }

  list_stats(angl)$prefiltered <- prefiltered_dt
  return(angl)
}

# ---------------------------------------------------------------------------
#' Select Genes Based on Statistical Thresholds from an anglemania_object
#'
#' Selects genes from an \code{\link{anglemania_object-class}} based on specified
#' thresholds for the absolute mean z-score and signal-to-noise ratio
#' (SNR) z-score. It updates the \code{integration_genes} slot of the
#' \code{anglemania_object} with the selected genes and associated
#' information.
#'
#' @param angl An \code{anglemania_object} containing statistical
#'   matrices such as mean z-scores and SNR z-scores.
#' @param zscore_mean_threshold Numeric value specifying the threshold for the
#'   absolute mean z-score. Default is 2.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   SNR z-score. Default is 2.
#' @param max_n_genes Integer specifying the maximum number of genes to select.
#'   If \code{NULL}, all genes that pass the thresholds are used. Default is
#'   \code{NULL}.
#' @return The input \code{anglemania_object} with the
#'   \code{integration_genes} slot updated to include the selected genes and
#'   their statistical information.
#' @importFrom stats quantile
#' @importFrom dplyr filter
#' @details
#' The function performs the following steps:
#' \enumerate{
#'  \item Checks if the input object is of class \code{\link{anglemania_object-class}}.
#'  \item If \code{max_n_genes} is not specified, it uses all genes that pass
#'     the thresholds.
#'  \item Identifies gene pairs where both the mean z-score and SNR z-score
#'     exceed the specified thresholds.
#'  \item If no gene pairs meet the criteria, it adjusts the thresholds to the
#'     99th percentile values of the corresponding statistics and re-selects.
#'  \item Extracts unique genes from the selected gene pairs using
#'     \code{\link{extract_rows_for_unique_genes}}.
#'  \item Updates the \code{integration_genes} slot of the
#'   \code{anglemania_object}
#'     with the selected genes and their statistics.
#' }
#' @examples
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl <- anglemania(angl)
#' angl <- select_genes(angl,
#'                       zscore_mean_threshold = 2.5,
#'                      zscore_sn_threshold = 2.5,
#'                      max_n_genes = 2000)
#' anglemania_genes <- get_anglemania_genes(angl)
#' # View the selected genes and use for integration
#' @seealso \code{\link{extract_rows_for_unique_genes}},
#'   \code{\link{intersect_genes}}, \code{\link{list_stats}}
#' @export
select_genes <- function(
    angl,
    zscore_mean_threshold = 2,
    zscore_sn_threshold = 2,
    max_n_genes = NULL) {
  if (!inherits(angl, "anglemania_object")) {
    stop("angl needs to be an anglemania_object")
  }

  if (is.null(max_n_genes)) {
    # If no max_n_genes specified, use all genes that pass the threshold
    max_n_genes <- length(intersect_genes(angl))
  }

  # Filter DF by thresholds
  filtered_genes_df <- list_stats(angl)$prefiltered %>%
    dplyr::filter(
      abs(mean_zscore) >= zscore_mean_threshold &
        sn_zscore >= zscore_sn_threshold
    )
  # Adjust thresholds if no genes passed the cutoff
  while (nrow(filtered_genes_df) == 0) {
    if (zscore_mean_threshold >= 0){
      message("No genes passed the cutoff.")
      zscore_mean_threshold <- zscore_mean_threshold - 0.1
        message(
          paste0(
            "Decreasing zscore_mean_threshold by 0.1: ",
            zscore_mean_threshold
          )
        )

      zscore_sn_threshold <- zscore_sn_threshold - 0.1
      message(
        paste0(
          "Decreasing zscore_sn_threshold by 0.1: ",
          zscore_sn_threshold
        )
      )
      filtered_genes_df <- list_stats(angl)$prefiltered %>%
        dplyr::filter(
          abs(mean_zscore) >= zscore_mean_threshold &
            sn_zscore >= zscore_sn_threshold
        )
    } else {
      stop("No genes passed the cutoff.")
    }
  }

  message(
    "If desired, you can re-run the selection of genes with a lower ",
    "zscore_mean_threshold and/or zscore_sn_threshold by using the ",
    "'select_genes' function. e.g.: angl <- select_genes(",
    "angl, zscore_mean_threshold = 1, zscore_sn_threshold = 1, ",
    "max_n_genes = 2000)"
  )
  message(
    "Please inspect get_anglemania_genes(angl)$info",
    " for info on the scores of the selected gene pairs."
  )



  # Order data frame
  filtered_genes_df <- filtered_genes_df[order(abs(filtered_genes_df$mean_zscore), decreasing = TRUE), ]
  angl@integration_genes$info <- filtered_genes_df
  selected_genes <- extract_rows_for_unique_genes(filtered_genes_df, max_n_genes)
  angl@integration_genes$genes <- 
    intersect_genes(angl)[selected_genes]
  print(paste0(
    "Selected ", length(selected_genes), " genes for integration."
  ))
  return(angl)
}


# ---------------------------------------------------------------------------
#' Adjusted big_cor function with different scaling function.
#' In comparison to big_scale(), we don't stop the scaling if there is a gene
#' with 0 variance.
#' #' Correlation
#'
#' Compute the correlation matrix of a Filebacked Big Matrix.
#' @param X An object of class [FBM][FBM-class].
#' @param backingfile The path to the backing file of the FBM.
#' @import bigstatsr
#' @return An FBM object with the correlation matrix of the data in X.
#' @examples
#' X <- bigstatsr::FBM(13, 17, init = rnorm(221))
#' angl_cor(X)
#' @export
angl_cor <- function(X,
                     backingfile = tempfile(tmpdir = getOption("FBM.dir"))) {
  # Adjusted scaling function: does not stop on 0 variance

  ind.row <- bigstatsr::rows_along(X)
  ind.col <- bigstatsr::cols_along(X)
  block.size <- bigstatsr::block_size(nrow(X))
  adj_big_scale <- function(center = TRUE, scale = TRUE) {
    function(X, ind.row, ind.col, ncores = 1) {
      bigstatsr:::check_args()
      m <- length(ind.col)

      if (center) {
        stats <- bigstatsr::big_colstats(X, ind.row, ind.col, ncores = ncores)
        means <- stats$sum / length(ind.row)
        sds <- if (scale) sqrt(stats$var) else rep(1, m)
        # Replace near-zero standard deviations with 1 (or a small constant)
        sds[sds < .Machine$double.eps] <- 1e-5
      } else {
        means <- rep(0, m)
        sds <- rep(1, m)
      }

      # Note: we do not check for near-zero sds here, so genes with zero variance are allowed.
      data.frame(center = means, scale = sds)
    }
  }

  # Scaling function for correlation that adjusts the standard deviation.
  cor.scaling <- function(X, ind.row, ind.col) {
    ms <- adj_big_scale(center = TRUE, scale = TRUE)(X, ind.row, ind.col)
    ms$scale <- ms$scale * sqrt(length(ind.row) - 1)
    ms
  }

  bigstatsr::big_crossprodSelf(X,
    fun.scaling = cor.scaling,
    ind.row = ind.row,
    ind.col = ind.col,
    block.size = block.size,
    backingfile = backingfile
  )
}
