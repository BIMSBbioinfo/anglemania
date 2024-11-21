# ---------------------------------------------------------------------------
# Compute Critical Angles Between Genes Across Samples in an anglemaniaObject
# ---------------------------------------------------------------------------

#' Compute Critical Angles Between Genes Across Samples in an anglemaniaObject
#'
#' @description
#' `anglemania` computes critical angles between genes across all samples
#' provided in an \code{\link{anglemaniaObject}}. It calculates angles,
#' transforms them to z-scores, computes statistical measures, and selects
#' genes based on specified thresholds.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Computes angles between genes for each sample in the
#'     \code{anglemania_object} using the specified \code{method}, via
#'     \code{\link{factorise}}.
#'   \item Transforms the angles to z-scores.
#'   \item Computes statistical measures (mean z-score, signal-to-noise ratio)
#'     across samples using \code{\link{get_list_stats}}.
#'   \item Selects genes based on specified z-score thresholds using
#'     \code{\link{select_genes}}.
#' }
#'
#' The computed statistics and selected genes are added to the
#' \code{anglemania_object}, which is returned.
#'
#' @param anglemania_object An \code{\link{anglemaniaObject}} containing gene
#' expression data and associated metadata.
#' @param method Character string specifying the method to use for calculating
#' the relationship between gene pairs. Default is \code{"pearson"}.
#' Other options include \code{"diem"}
#' (see \url{https://bytez.com/docs/arxiv/2407.08623/paper}).
#' @param zscore_mean_threshold Numeric value specifying the threshold
#'  for the mean z-score. Default is \code{2.5}.
#' @param zscore_sn_threshold Numeric value specifying the threshold for the
#'   signal-to-noise z-score. Default is \code{2.5}.
#' @param max_n_genes Integer specifying the maximum number of genes to select.
#'   Default is \code{2000}.
#'
#' @return An updated \code{\link{anglemaniaObject}} with computed statistics and
#'   selected genes based on the specified thresholds.
#'
#' @importFrom pbapply pblapply
#' @importFrom pbapply pboptions
#'
#' @seealso
#'   \code{\link{create_anglemaniaObject}},
#'   \code{\link{get_list_stats}},
#'   \code{\link{select_genes}},
#'   \code{\link{factorise}},
#'   \code{\link[bigstatsr]{big_apply}},
#'   \url{https://arxiv.org/abs/1306.0256}
#'
#' @examples
#' \donttest{
#' load(system.file(
#'  "extdata",
#'  "seurat_splatter_sim.RData",
#'  package = "anglemania"))
#' 
#' angl <- create_anglemaniaObject(se,
#'  batch_key = batch_key,
#'  min_cells_per_gene = 1
#'  )
#'
#' angl <- anglemania(
#'   angl,
#'   method = "pearson",
#'   zscore_mean_threshold = 2,
#'   zscore_sn_threshold = 2,
#'   max_n_genes = 2000
#' )
#'
#' # Access the selected genes
#' selected_genes <- get_anglemania_genes(angl)
#' }
#' selected_genes[1:10]
#' @export
anglemania <- function(
    anglemania_object,
    method = "pearson",
    zscore_mean_threshold = 2.5,
    zscore_sn_threshold = 2.5,
    max_n_genes = 2000) {
  # Validate inputs
  if (!inherits(anglemania_object, "anglemaniaObject")) {
    stop("anglemania_object needs to be an anglemaniaObject")
  }

  # dataset_key and batch_key are checked in create_anglemaniaObject!

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
  matrix_list(anglemania_object) <- pbapply::pblapply(
    matrix_list(anglemania_object),
    function(x) {
      factorise(
        x_mat = x,
        method = method,
        seed = 1
      )
    }
  )

  message("Computing statistics...")
  list_stats(anglemania_object) <- get_list_stats(anglemania_object)
  invisible(gc())

  message("Filtering features...")
  anglemania_object <- select_genes(
    anglemania_object,
    zscore_mean_threshold = zscore_mean_threshold,
    zscore_sn_threshold = zscore_sn_threshold,
    max_n_genes = max_n_genes
  )

  return(anglemania_object)
}
