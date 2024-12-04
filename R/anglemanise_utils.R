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
#' This function takes an \code{anglemaniaObject} containing a list of
#' \code{\link[bigstatsr]{FBM}}s and calculates the mean of every element.
#' If the list is empty or the FBMs have different dimensions,
#' it throws an error.
#'
#' @param anglemania_object An \code{anglemaniaObject} containing the list of FBMs.
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
#' # Construct the anglemaniaObject
#' anglemania_object <- new(
#'   "anglemaniaObject",
#'   weights = weights,
#'   matrix_list = fbm_list
#' )
#' big_mat_list_mean(anglemania_object)
#' @export
big_mat_list_mean <- function(anglemania_object) {
  if (!inherits(anglemania_object, "anglemaniaObject")) {
    stop("anglemania_object needs to be an anglemaniaObject")
  }

  # Check if all matrices in the list have the same dimensions
  if (!all(sapply(
    matrix_list(anglemania_object),
    function(x) {
      identical(dim(matrix_list(anglemania_object)[[1]]), dim(x))
    }
  ))) {
    stop("All matrices in the list must have the same dimensions.")
  }

  n_col <- ncol(matrix_list(anglemania_object)[[1]])
  n_row <- nrow(matrix_list(anglemania_object)[[1]])
  mat_mean_zscore <- bigstatsr::FBM(n_row, n_col)

  bigstatsr::big_apply(
    mat_mean_zscore,
    a.FUN = function(X, ind) {
      X.sub <- X[, ind, drop = FALSE]
      wrap_mean <- function(final_mat, batch) {
        batch_mat <- matrix_list(anglemania_object)[[batch]][
          , ind,
          drop = FALSE
        ] * anglemania_object@weights[batch]
        final_mat <- final_mat + batch_mat
      }
      # Run function in Reduce statement on names of weights vector
      X.sub <- Reduce(
        wrap_mean,
        names(anglemania_object@weights),
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
#' for each element across a list of FBMs in an \code{anglemaniaObject}.
#'
#' @param anglemania_object An \code{anglemaniaObject} containing the list of FBMs.
#' @return A list containing three matrices: \code{mean_zscore},
#'   \code{sds_zscore}, and \code{sn_zscore}.
#' @importFrom bigstatsr FBM big_apply
#' @examples
#' load(system.file(
#' "extdata",
#'  "seurat_splatter_sim.RData",
#'  package = "anglemania")
#' )
#' anglemania_object <- create_anglemaniaObject(se, batch_key = "Batch")
#' anglemania_object <- anglemania(anglemania_object)
#' list_stats(anglemania_object) <- get_list_stats(anglemania_object)
#' str(list_stats(anglemania_object))
#' @seealso \code{\link[bigstatsr]{big_apply}}, \code{\link[bigstatsr]{FBM}}
#' @export
get_list_stats <- function(anglemania_object) {
  if (!inherits(anglemania_object, "anglemaniaObject")) {
    stop("anglemania_object needs to be an anglemaniaObject")
  }

  # Check if all matrices in the list have the same dimensions
  if (!all(sapply(
    matrix_list(anglemania_object),
    function(x) {
      identical(dim(matrix_list(anglemania_object)[[1]]), dim(x))
    }
  ))) {
    stop("All matrices in the list need to have the same dimensions.")
  }

  message("Weighting matrix_list...")
  message("Calculating mean...")

  mat_mean_zscore <- big_mat_list_mean(anglemania_object)

  n_col <- ncol(matrix_list(anglemania_object)[[1]])
  n_row <- nrow(matrix_list(anglemania_object)[[1]])
  mat_sds_zscore <- bigstatsr::FBM(n_row, n_col)

  message("Calculating sds...")

  bigstatsr::big_apply(
    mat_sds_zscore,
    a.FUN = function(X, ind) {
      wrap_sds <- function(final_mat, batch) {
        batch_mat <- matrix_list(anglemania_object)[[batch]]
        final_mat <- final_mat + (
          batch_mat[, ind, drop = FALSE] -
            mat_mean_zscore[, ind, drop = FALSE]
        )^2 * anglemania_object@weights[batch]
      }
      X.sub <- X[, ind, drop = FALSE]
      X.sub <- Reduce(
        wrap_sds,
        names(anglemania_object@weights),
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
    mean_zscore = mat_mean_zscore[],
    sds_zscore  = mat_sds_zscore[],
    sn_zscore   = mat_sn_zscore[]
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
#' Select Genes Based on Statistical Thresholds from an anglemaniaObject
#'
#' Selects genes from an \code{\link{anglemaniaObject-class}} based on specified
#' thresholds for the absolute mean z-score and signal-to-noise ratio
#' (SNR) z-score. It updates the \code{integration_genes} slot of the
#' \code{anglemaniaObject} with the selected genes and associated
#' information.
#'
#' @param anglemania_object An \code{anglemaniaObject} containing statistical
#'   matrices such as mean z-scores and SNR z-scores.
#' @param zscore_mean_threshold Numeric value specifying the threshold for the
#'   absolute mean z-score. Default is 2.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   SNR z-score. Default is 2.
#' @param max_n_genes Integer specifying the maximum number of genes to select.
#'   If \code{NULL}, all genes that pass the thresholds are used. Default is
#'   \code{NULL}.
#' @param direction whether to select genes with positive, negative or both mean z-scores. Default is "both"
#' @param adjust_threshols whether to automatically adjust threholds if the selected genes do not meet the thresholds. Default is TRUE

#' @return The input \code{anglemaniaObject} with the
#'   \code{integration_genes} slot updated to include the selected genes and
#'   their statistical information.
#' @importFrom stats quantile
#' @details
#' The function performs the following steps:
#' \enumerate{
#'  \item Checks if the input object is of class \code{\link{anglemaniaObject-class}}.
#'  \item If \code{max_n_genes} is not specified, it uses all genes that pass
#'     the thresholds.
#'  \item Identifies gene pairs where both the mean z-score and SNR z-score
#'     exceed the specified thresholds.
#'  \item If no gene pairs meet the criteria, it adjusts the thresholds to the
#'     99th percentile values of the corresponding statistics and re-selects.
#'  \item Extracts unique genes from the selected gene pairs using
#'     \code{\link{extract_rows_for_unique_genes}}.
#'  \item Updates the \code{integration_genes} slot of the
#'   \code{anglemaniaObject}
#'     with the selected genes and their statistics.
#' }
#' @examples
#' load(system.file(
#'  "extdata",
#'  "seurat_splatter_sim.RData",
#'  package = "anglemania"))
#' 
#' anglemania_object <- create_anglemaniaObject(se,
#'  batch_key = batch_key,
#'  min_cells_per_gene = 1
#'  )
#'
#' anglemania_object <- anglemania(
#'   anglemania_object,
#'   method = "pearson",
#'   zscore_mean_threshold = 2,
#'   zscore_sn_threshold = 2,
#'   max_n_genes = 2000
#' )
#' anglemania_object <- select_genes(anglemania_object,
#'                       zscore_mean_threshold = 2.5,
#'                      zscore_sn_threshold = 2.5,
#'                      max_n_genes = 2000)
#' anglemania_genes <- get_anglemania_genes(anglemania_object)
#' # View the selected genes and use for integration
#' @seealso \code{\link{extract_rows_for_unique_genes}},
#'   \code{\link{intersect_genes}}, \code{\link{list_stats}}
#' @export
select_genes <- function(
    anglemania_object,
    zscore_mean_threshold = 2,
    zscore_sn_threshold = 2,
    max_n_genes = NULL,
    direction   = "both",
    adjust_thresholds = TRUE
) {
  if (!inherits(anglemania_object, "anglemaniaObject")) {
    stop("anglemania_object needs to be an anglemaniaObject")
  }
  if (!direction %in% c("both","positive","negative"))
    stop("direction can only be 'both', 'positive' or 'negative'")
  

  if (is.null(max_n_genes)) {
    # If no max_n_genes specified, use all genes that pass the threshold
    max_n_genes <- length(intersect_genes(anglemania_object))
  }

  # Selects the direction of conserved genes
  if(direction == "both"){
    gene_ind <- which(
      upper.tri(list_stats(anglemania_object)$sn_zscore) &
        (list_stats(anglemania_object)$sn_zscore >= zscore_sn_threshold) &
        (abs(list_stats(anglemania_object)$mean_zscore) >= zscore_mean_threshold),
      arr.ind = TRUE
    )
  }else if(direction == "positive"){
    gene_ind <- which(
      upper.tri(list_stats(anglemania_object)$sn_zscore) &
        (list_stats(anglemania_object)$sn_zscore >= zscore_sn_threshold) &
        (list_stats(anglemania_object)$mean_zscore >= zscore_mean_threshold),
      arr.ind = TRUE
    )
  }else if(direction == "negative"){
    gene_ind <- which(
      upper.tri(list_stats(anglemania_object)$sn_zscore) &
        (list_stats(anglemania_object)$sn_zscore >= zscore_sn_threshold) &
        (list_stats(anglemania_object)$mean_zscore <= zscore_mean_threshold),
      arr.ind = TRUE
    )
  }

  # Adjust thresholds if no genes passed the cutoff
  if (nrow(gene_ind) == 0 && adjust_thresholds) {
    message("No genes passed the cutoff.")
    quantile95mean <- stats::quantile(
      abs(list_stats(anglemania_object)$mean_zscore),
      0.95,
      na.rm = TRUE
    )
    quantile95sn <- stats::quantile(
      list_stats(anglemania_object)$sn_zscore,
      0.95,
      na.rm = TRUE
    )
    if (quantile95mean < zscore_mean_threshold) {
      zscore_mean_threshold <- quantile95mean
      message(
        paste0(
          "zscore_mean_threshold is lower than the 95% quantile of the ",
          "absolute mean z-scores. Setting zscore_mean_threshold to: ",
          zscore_mean_threshold
        )
      )
    }
    if (quantile95sn < zscore_sn_threshold) {
      zscore_sn_threshold <- quantile95sn
      message(
        paste0(
          "zscore_sn_threshold is higher than 95% quantile. Setting ",
          "zscore_sn_threshold to: ",
          round(zscore_sn_threshold, 2)
        )
      )
    }
    print(
      paste0(
        "zscore_mean_threshold: ", zscore_mean_threshold,
        " zscore_sn_threshold: ", round(zscore_sn_threshold, 2)
      )
    )

    # Re-run the selection
    gene_ind <- which(
      upper.tri(list_stats(anglemania_object)$sn_zscore) &
        (list_stats(anglemania_object)$sn_zscore >= zscore_sn_threshold) &
        (abs(list_stats(anglemania_object)$mean_zscore) >= zscore_mean_threshold),
      arr.ind = TRUE
    )

    message(
      "If desired, you can re-run the selection of genes with a lower ",
      "zscore_mean_threshold and/or zscore_sn_threshold by using the ",
      "'select_genes' function. e.g.: anglemania_object <- select_genes(",
      "anglemania_object, zscore_mean_threshold = 1, zscore_sn_threshold = 1, ",
      "max_n_genes = 2000)"
    )
    message(
      "Please inspect get_anglemania_genes(anglemania_object)$info",
      " for info on the scores of the selected gene pairs."
    )
  }

  top_n = data.frame(
    geneA = gene_ind[, 1],
    geneB = gene_ind[, 2],
    zscore  = list_stats(anglemania_object)$mean_zscore[gene_ind],
    snscore = list_stats(anglemania_object)$sn_zscore[gene_ind],
    gene_nameA = intersect_genes(anglemania_object)[gene_ind[, 1]],
    gene_nameB = intersect_genes(anglemania_object)[gene_ind[, 2]]
  )

  # Order data frame
  top_n <- top_n[order(abs(top_n$zscore), decreasing = TRUE), ]
  anglemania_object@integration_genes$info <- top_n
  selected_genes <- extract_rows_for_unique_genes(top_n, max_n_genes)
  anglemania_object@integration_genes$genes <- 
    intersect_genes(anglemania_object)[selected_genes]

  return(anglemania_object)
}
