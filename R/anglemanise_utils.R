# ---------------------------------------------------------------------------
# Utility functions for the anglemania package
# ---------------------------------------------------------------------------
#' @title Utility Functions for the anglemania Package
#' @description A collection of utility functions used within the 
#' \pkg{anglemania} package
#' for manipulating FBMs, calculating statistics, and selecting genes.
#' @name anglemanise_utils
#' @rdname anglemanise_utils
#' @keywords internal
NULL

#' @describeIn anglemanise_utils Convert a sparse matrix into a file-backed matrix
#' (\code{\link[bigstatsr]{FBM}}) with efficient memory usage.
#'
#' @param s_mat A sparse matrix.
#' @return An \code{\link[bigstatsr]{FBM}} object from the \pkg{bigstatsr}
#'   package.
#' @importFrom bigstatsr FBM
#' @examples
#' s_mat <- Matrix::rsparsematrix(nrow = 10, ncol = 5, density = 0.3)
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
#' @describeIn anglemanise_utils Compute mean, standard deviation, variance,
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
      stat_value = replace_with_na(stat_value)
      stat_value
    })
  })

  # Calculate standard deviation
  dstat$sd = sqrt(dstat$var)

  return(dstat)
}

# ---------------------------------------------------------------------------
#' @describeIn anglemanise_utils Calculates the element-wise mean from a list
#' of \code{\link[bigstatsr]{FBM}}s stored in an \code{anglemania_object}.
#'
#' @param angl An \code{anglemania_object} containing the list of FBMs.
#'   In this case, the FBMs are the angle matrices computed in \code{factorise}.
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
  if (
    !all(sapply(
      matrix_list(angl),
      function(x) {
        identical(dim(matrix_list(angl)[[1]]), dim(x))
      }
    ))
  ) {
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
        batch_mat <- matrix_list(angl)[[batch]][,
          ind,
          drop = FALSE
        ]

        # calculates the number of samples in which the feature is present, weighted
        # by the dataset weight
        nmat <- (!is.na(batch_mat) + 0) * angl@weights[batch]
        final_mat <- final_mat + batch_mat
        final_mat[is.na(final_mat)] = 0
        list(final_mat, nmat)
      }
      # Run function in Reduce statement on names of weights vector
      lmats <- lapply(names(angl@weights), function(batch) {
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
#' @describeIn anglemanise_utils Calculate mean, standard deviation, and SNR
#' across a list of FBMs stored in an \code{anglemania_object}.
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
  if (
    !all(sapply(
      matrix_list(angl),
      function(x) {
        identical(dim(matrix_list(angl)[[1]]), dim(x))
      }
    ))
  ) {
    stop("All matrices in the list need to have the same dimensions.")
  }

  message("Weighting matrix_list...")
  message("Calculating mean...")
  mat_mean_zscore <- big_mat_list_mean(angl)

  n_col <- ncol(matrix_list(angl)[[1]])
  n_row <- nrow(matrix_list(angl)[[1]])
  mat_sds_zscore <- bigstatsr::FBM(n_row, n_col)

  message("Calculating sds...")
  # calculates the weighted sds for the big matrix
  bigstatsr::big_apply(
    mat_sds_zscore,
    a.FUN = function(X, ind) {
      X.sub <- X[, ind, drop = FALSE]
      wrap_mean <- function(final_mat, batch) {
        batch_mat <- matrix_list(angl)[[batch]][,
          ind,
          drop = FALSE
        ]

        # calculates the number of samples in which the feature is present, weighted
        # by the dataset weight
        nmat <- (!is.na(batch_mat) + 0) * angl@weights[batch]

        # calculates the weighted standard deviation
        final_mat <- final_mat +
          (batch_mat - mat_mean_zscore[, ind, drop = FALSE])^2 *
            angl@weights[batch]
        final_mat[is.na(final_mat)] = 0
        list(final_mat, nmat)
      }
      # Run function in Reduce statement on names of weights vector
      lmats <- lapply(names(angl@weights), function(batch) {
        wrap_mean(X.sub, batch)
      })
      # this sums up the z-scores across samples
      m_sd <- Reduce("+", lapply(lmats, '[[', 1))
      # this gets the number of samples in which the feature was present
      w_m_sum    <- Reduce("+", lapply(lmats, function(x)x[[2]]  ))
      w_m_sum_sq <- Reduce("+", lapply(lmats, function(x)x[[2]]^2))
      
      # creates an unbiased esitmator for a sample based wighted variance calculation
      denominator = w_m_sum - (w_m_sum_sq / w_m_sum)

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

# ---------------------------------------------------------------------------
#' @describeIn anglemanise_utils Extract unique gene identifiers from gene pairs,
#' returning up to a specified maximum number.
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
#' @describeIn anglemanise_utils Preselect genes from an 
#' \code{anglemania_object}
#' based on z-score thresholds, to simplify downstream filtering.
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
    zscore_sn_threshold = 1
  ) {
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
#' @describeIn anglemanise_utils Select genes from an \code{anglemania_object}
#' based on z-score thresholds for mean and signal-to-noise ratio (SNR).
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
#' @param direction whether to select genes with positive, negative or both
#' mean z-scores. Default is "both"
#' @param adjust_thresholds whether to automatically adjust threholds if the
#' selected genes do not meet the thresholds. Default is TRUE
#' @return The input \code{anglemaniaObject} with the
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
  max_n_genes = NULL,
  direction = "both",
  adjust_thresholds = TRUE
) {
  if (!inherits(angl, "anglemania_object")) {
    stop("anglemania_object needs to be an anglemania_object")
  }
  if (!direction %in% c("both", "positive", "negative"))
    stop("direction can only be 'both', 'positive' or 'negative'")

  if (is.null(max_n_genes)) {
    # If no max_n_genes specified, use all genes that pass the threshold
    max_n_genes <- length(intersect_genes(angl))
  }

  # Selects the direction of conserved genes
  prefiltered = list_stats(angl)$prefiltered
  if (direction == "both") {
    filtered_genes_df = subset(
      prefiltered,
      sn_zscore >= zscore_sn_threshold &
        abs(prefiltered$mean_zscore) >= zscore_mean_threshold
    )
  } else if (direction == "positive") {
    filtered_genes_df = subset(
      prefiltered,
      sn_zscore >= zscore_sn_threshold &
        prefiltered$mean_zscore >= zscore_mean_threshold
    )
  } else if (direction == "negative") {
    filtered_genes_df = subset(
      prefiltered,
      sn_zscore >= zscore_sn_threshold &
        prefiltered$mean_zscore <= zscore_mean_threshold
    )
  }

  # Adjust thresholds if no genes passed the cutoff
  while (nrow(filtered_genes_df) == 0 && adjust_thresholds) {
    message("No genes passed the cutoff.")
    zscore_sn_threshold <- zscore_sn_threshold - 0.1
    zscore_mean_threshold <- zscore_mean_threshold - 0.1
    message(
      paste0(
        "Decreasing zscore mean and sn thresholds by 0.1: \n",
        "zscore_mean_threshold: ", round(zscore_mean_threshold, 2), "\n",
        "zscore_sn_threshold: ", round(zscore_sn_threshold, 2)
      )
    )
    filtered_genes_df <- list_stats(angl)$prefiltered %>%
      dplyr::filter(
        abs(mean_zscore) >= zscore_mean_threshold &
            sn_zscore >= zscore_sn_threshold
        )
  }

  filtered_genes_df$gene_nameA = intersect_genes(angl)[filtered_genes_df[, 1]]
  filtered_genes_df$gene_nameB = intersect_genes(angl)[filtered_genes_df[, 2]]
  # message(
  #   "Please inspect get_anglemania_genes(angl)$info",
  #   " for info on the scores of the selected gene pairs."
  # )

  # Order data frame
  filtered_genes_df <- filtered_genes_df[
    order(abs(filtered_genes_df$mean_zscore), decreasing = TRUE),
  ]
  angl@integration_genes$info <- filtered_genes_df
  selected_genes <- extract_rows_for_unique_genes(
    filtered_genes_df,
    max_n_genes
  )
  angl@integration_genes$genes <-
    intersect_genes(angl)[selected_genes]
  print(paste0(
    "Selected ",
    length(selected_genes),
    " genes for integration."
  ))
  return(angl)
}

# ---------------------------------------------------------------------------
#' @describeIn anglemanise_utils replace Nan and Inf values with NA
#'
#' @description
#' Replaces all NaN and Inf values in a numeric vector with NA.
#' @param v A numeric vector.
#' @return A numeric vector with NaN and Inf values replaced with NA.
#' @examples
#' v <- c(1, 2, 3, 4, 5, 6, 7, Inf, 9, NA)
#' v <- replace_with_na(v)
#' v
#' @export
replace_with_na = function(v) {
  v[is.nan(v)] = NA
  v[is.infinite(v)] = NA
  v
}



# ---------------------------------------------------------------------------
#' @describeIn anglemanise_utils Overlap of Variable Genes Across Thresholds.
#' This function computes the overlap between variable genes identified in a
#' Seurat object and the gene sets selected from an "anglemania" analysis
#' object (`ang_obj`) at multiple threshold levels. For each threshold, genes
#' are selected from `ang_obj` by applying `select_genes`, and the intersection
#' of these genes with the variable features of the Seurat object is calculated.
#' The result is a summary data frame that includes the number and percentage of
#' variable genes overlapping with the anglemania-selected genes, as well as the
#' number of genes selected at each threshold.
#' @param seu A Seurat object. If variable features have not been identified,
#'   they will be computed using \code{\link[Seurat]{FindVariableFeatures}}.
#' @param angl An anglemania analysis object. This is used by `select_genes`
#'  and `get_anglemania_genes` to select genes based on the provided thresholds.
#' @param zscore_mean_thresholds A numeric vector of zscore mean thresholds used
#'   for gene selection from `ang_obj`. Default is \code{c(0.5, 1:15)}.
#' @param zscore_sn_thresholds A numeric vector of zscore signal to noise ratio
#'   thresholds used for gene selection from `ang_obj`. Default is
#'   \code{c(0.5, 1:15)}.
#' @param adjust_thresholds Logical. Whether to use adaptive thresholding in
#'   `select_genes`. Default is FALSE.
#' @param layer Character. The assay layer in the Seurat object to use when
#'   calculating variable features (e.g., "counts", "data", or "logcounts").
#'   Default is \code{"data"}.
#' @importFrom Seurat FindVariableFeatures
#' @importFrom SeuratObject VariableFeatures
#' @return A data frame with the following columns:
#' \describe{
#'   \item{threshold}{The threshold value used for gene selection.}
#'   \item{intersecting_genes}{The number of genes that are both variable in
#'     \code{seu} and selected by \code{ang_obj}.}
#'   \item{perc_var_genes}{The proportion of variable genes in \code{seu} that
#'     are also selected by \code{ang_obj}.}
#'   \item{number_angl_genes}{The total number of genes selected by
#'     \code{ang_obj} at the given threshold.}
#' }
#' @examples
#' sce <- sce_example()
#' seu <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(sce))
#' seu <- Seurat::FindVariableFeatures(seu)
#' angl <- create_anglemania_object(
#'   sce,
#'   batch_key = "batch",
#'   min_cells_per_gene = 1
#' )
#' angl <- anglemania(angl)
#' head(variable_genes_overlap(seu, angl))
#' @export
variable_genes_overlap = function(
  seu,
  angl,
  zscore_mean_thresholds = NULL,
  zscore_sn_thresholds = NULL,
  adjust_thresholds = FALSE,
  layer = "data",
  length_out = 6
) {
  # Check if the Seurat object has variable features. If not, compute them.
  if (length(SeuratObject::VariableFeatures(seu)) == 0) {
    seu_list <- Seurat::SplitObject(
      seu,
      split.by = batch_key(angl)
    )
    hvgs <- lapply(seu_list, function(x) Seurat::NormalizeData(x)) %>%
      Seurat::SelectIntegrationFeatures()
  }

  if (is.null(zscore_mean_thresholds)) {
    zscore_mean_thresholds <- seq(
      0.5,
      max(list_stats(angl)$mean_zscore[], na.rm = TRUE),
      length.out = length_out
    )
  }
  if (is.null(zscore_sn_thresholds)) {
    zscore_sn_thresholds <- seq(
      0.5,
      max(list_stats(angl)$sn_zscore[], na.rm = TRUE),
      length.out = length_out
    )
  }

  # For each threshold, select genes from anglemania object and extract their names.
  dparams <- expand.grid(zscore_mean_thresholds, zscore_sn_thresholds)
  colnames(dparams) <- c("zscore_mean_thresholds", "zscore_sn_thresholds")
  lg <- lapply(seq_len(nrow(dparams)), function(i) {
    angl <- select_genes(
      angl,
      zscore_mean_threshold = dparams$zscore_mean_thresholds[i],
      zscore_sn_threshold = dparams$zscore_sn_thresholds[i],
      direction = "both",
      adjust_thresholds = adjust_thresholds
    )
    gg <- get_anglemania_genes(angl)
    if (length(gg) < 0) gg <- ""
    return(gg)
  })

  # Compute the intersection of anglemania-selected genes with variable
  # features of the Seurat object.
  lg_int <- lapply(
    lg,
    function(x) length(intersect(x, hvgs))
  )

  # Construct a data frame summarizing the intersections and percentages.
  dg_int <- dparams %>%
    data.frame() %>%
    dplyr::mutate(intersecting_genes = unlist(lg_int)) %>%
    dplyr::mutate(
      perc_var_genes = intersecting_genes /
        length(hvgs)
    ) %>%
    dplyr::mutate(number_angl_genes = sapply(lg, length)) %>%
    dplyr::mutate(perc_ang_genes = intersecting_genes / number_angl_genes)

  return(dg_int)
}

# ---------------------------------------------------------------------------
#' @describeIn anglemanise_utils normalize matrix
#' Normalize a Filebacked Big Matrix (`FBM`) using either total-count
#' scaling or residuals from a linear model. Intended for use with
#' single-cell RNA-seq gene expression data.
#'
#' @param x_mat A `bigstatsr::FBM` object containing the matrix to normalize
#'   (typically genes x cells).
#' @param normalization_method A character string specifying the normalization
#'   method to use. One of `"divide_by_total_counts"` (default) or
#'   `"find_residuals"`. 
#'   - `"divide_by_total_counts"` normalizes each cell by its total expression
#'     count and applies log1p.
#'   - `"find_residuals"` computes log1p-transformed residuals after regressing
#'     out total expression.
#'
#' @return The input `FBM` object with normalized values written back in place. 
#'   This function modifies the input `x_mat` by reference.
#'
#' @importFrom bigstatsr big_apply
#' @importFrom checkmate assertChoice assertClass
#' @examples
#' library(bigstatsr)
#' set.seed(42)
#' mat <- matrix(rpois(1000, lambda = 5), nrow = 100, ncol = 10)
#' fbm <- as_FBM(mat)
#'
#' normalize_matrix(fbm, normalization_method = "divide_by_total_counts")[1:5, 1:5]
#' normalize_matrix(fbm, normalization_method = "find_residuals")[1:5, 1:5]
#'
#' @export
normalize_matrix = function(
  x_mat,
  normalization_method = "divide_by_total_counts"
){
  checkmate::assertChoice(normalization_method, c("divide_by_total_counts", "find_residuals"))
  checkmate::assertClass(x_mat, "FBM")

  bigstatsr::big_apply(x_mat, a.FUN = function(X, ind) {
    X.sub <- X[, ind, drop = FALSE]
    if(normalization_method == "divide_by_total_counts"){
      # Normalize the data:
      #   divide gene counts by the number of total counts per cell,
      #   multiply by 10,000 (scaling factor like in Seurat)
      # +1 prevents NaN values when working with partial features
    
      X.sub <- t(t(X.sub) / (colSums(X.sub)) * 10000)
      X.res <- log1p(X.sub)
    } else if(normalization_method == "find_residuals"){
      total <- log1p(colSums(X.sub))
      X.sub <- log1p(X.sub)
      X.res = t(apply(X.sub, 1, function(x)residuals(lm(x ~ total))))
    }

    X[, ind] <- X.res
    NULL
  })
  x_mat
}