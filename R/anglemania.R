# ---------------------------------------------------------------------------
# Compute Critical Angles Between Genes Across Samples in an anglemania_object
# ---------------------------------------------------------------------------

#' Compute Critical Angles Between Genes Across Samples in an anglemania_object
#'
#' @description
#' `anglemania` computes critical angles between genes across all samples
#' provided in an \code{\link{anglemania_object-class}}. It calculates angles,
#' transforms them to z-scores, computes statistical measures, and selects
#' genes based on specified thresholds.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Computes angles between genes for each sample in the
#'     \code{angl} using the specified \code{method}, via
#'     \code{\link{factorise}}.
#'   \item Transforms the angles to z-scores.
#'   \item Computes statistical measures (mean z-score, signal-to-noise ratio)
#'     across samples using \code{\link{get_list_stats}}.
#'   \item Selects genes based on specified z-score thresholds using
#'     \code{\link{select_genes}}.
#' }
#'
#' The computed statistics and selected genes are added to the
#' \code{angl}, which is returned.
#'
#' @param angl An \code{\link{anglemania_object-class}} containing gene
#' expression data and associated metadata.
#' @param method Character string specifying the method to use for calculating
#' the relationship between gene pairs. Default is \code{"cosine"}.
#' Other options include \code{"diem"}
#' (see \url{https://bytez.com/docs/arxiv/2407.08623/paper}).
#' @param zscore_mean_threshold Numeric value specifying the threshold
#'  for the mean z-score. Default is \code{2.5}.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   signal-to-noise z-score. Default is \code{2.5}.
#' @param max_n_genes Integer specifying the maximum number of genes to select.
#' @param permute_row_or_column Character "row" or "column", whether permutations should be executed row-wise or column wise. Default is \code{"column"}
#' @param permutation_function Character "sample" or "permute_nonzero". If sample, then sample is used for constructing background distributions. If permute_nonzero, then only non-zero values are permuted. Default is \code{"sample"}
#'   Default is \code{2000}.
#'
#' @return An updated \code{\link{anglemania_object-class}} with computed statistics and
#'   selected genes based on the specified thresholds.
#'
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
#'
#' @seealso
#'   \code{\link{create_anglemania_object}},
#'   \code{\link{get_list_stats}},
#'   \code{\link{select_genes}},
#'   \code{\link{factorise}},
#'   \code{\link[bigstatsr]{big_apply}},
#'   \url{https://arxiv.org/abs/1306.0256}
#'
#' @examples
#' # Set seed for reproducibility (optional)
#' sce <- sce_example()
#' angl <- create_anglemania_object(sce, batch_key = "batch")
#' angl <- anglemania(
#'   angl,
#'   method = "cosine",
#'   zscore_mean_threshold = 2,
#'   zscore_sn_threshold = 2,
#'   max_n_genes = 2000
#' )
#'
#' # Access the selected genes
#' selected_genes <- get_anglemania_genes(angl)
#' selected_genes[1:10]
#' @export
anglemania <- function(
    angl,
    zscore_mean_threshold = 2.5,
    zscore_sn_threshold = 2.5,
    max_n_genes = NULL,
    method = "cosine",
    permute_row_or_column = "columns",
    permutation_function = "sample"
) {
  # Validate inputs
  if (!inherits(angl, "anglemania_object")) {
    stop("angl needs to be an anglemania_object")
  }

  # dataset_key and batch_key are checked in create_anglemania_object!

  if (!is.numeric(zscore_mean_threshold) || zscore_mean_threshold < 0) {
    stop("zscore_mean_threshold has to be a non-negative number")
  }
  if (!is.numeric(zscore_sn_threshold) || zscore_sn_threshold < 0) {
    stop("zscore_sn_threshold has to be a non-negative number")
  }
  if (!is.null(max_n_genes) &&
    (!is.numeric(max_n_genes) || max_n_genes < 1)) {
    stop("max_n_genes has to be a positive integer")
  }

  # Process inputs
  pbapply::pboptions(
    type = "timer",
    style = 1,
    char = "=",
    title = "anglemania"
  )
  message("Computing angles and transforming to z-scores...")
  matrix_list(angl) <- pbapply::pblapply(
    matrix_list(angl),
    function(x) {
      factorise(
        x_mat = x,
        method = method,
        seed = 1,
        permute_row_or_column = permute_row_or_column,
        permutation_function  = permutation_function
      )
    }
  )

  message("Computing statistics...")
  list_stats(angl) <- get_list_stats(angl)
  invisible(gc())

  # this corrects sn values when sds is zero - can happen if only one pair of the gene is found 
  # in one sample
  list_stats(anglemania_object)$sn_zscore[list_stats(anglemania_object)$sds_zscore == 0] = 0

  message("Filtering features...")
  angl <- prefilter_angl(
    angl,
    zscore_mean_threshold = 1,
    zscore_sn_threshold = 1
  )

  message("Selecting features...")
  angl <- select_genes(
    angl,
    zscore_mean_threshold = zscore_mean_threshold,
    zscore_sn_threshold = zscore_sn_threshold,
    max_n_genes = max_n_genes
  )

  return(angl)
}
