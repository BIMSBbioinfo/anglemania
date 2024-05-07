# --------------------------------------------------------------------------------------------- #
# Utility functions for the anglemanise package
# --------------------------------------------------------------------------------------------- #
#' Convert a sparse matrix into an FBM
#'
#' Converts a sparse matrix into an FBM with efficient memory usage.
#'
#' @param s_mat A sparse matrix.
#' @return A FBM object from the bigstatsr package.
#' @importFrom bigstatsr FBM
sparse_to_fbm <- function(s_mat) {
    n <- nrow(s_mat)
    p <- ncol(s_mat)
    X <- bigstatsr::FBM(n, p)
    big_apply(X, a.FUN = function(X, ind) {
        X[, ind] <- s_mat[, ind] %>% as.matrix()
        NULL
    }, a.combine = "c", block.size = 200)

    return(X)
}

# --------------------------------------------------------------------------------------------- #
#' Permute counts in a Seurat object
#'
#' This function permutes the counts in the assay data of a Seurat object either by columns or rows.
#' @param seu A Seurat object.
#' @param which A character string specifying whether to permute "columns" or "rows".
#' @return A new Seurat object with permuted count data.
#' @examples
#' seu_permuted = permute_counts(seu, which = "columns")
permute_counts <- function(X.sub, which = "columns", seed=1) {
    set.seed(seed)
    if (which == "columns") {
        X.sub <- apply(X.sub, 2, sample)
    } else if (which == "rows") {
        X.sub <- t(apply(X.sub, 1, sample))
    }
}

# --------------------------------------------------------------------------------------------- #
#' Compute mean and standard deviation of the correlation matrix
#'
#' Computes the mean and standard deviation of the correlation matrix
#' using the big_apply function.
#'
#' @param corr_matrix A FBM object from the bigstatsr package.
#' @return A list with two entries: \code{mean} and \code{sd}.
#' @importFrom bigstatsr big_apply
get_dstat <- function(corr_matrix) {
    # check if FBM is used as input
    if (!inherits(corr_matrix, "FBM")) {
        stop("corr_matrix has to be a FBM object")
    }
    # calculate mean and standard deviation
    # n = number of entries that are not NA ==> so that when calling mean or sd, we do not count NA values, 
    #   corresponding to mean(..., na.rm = TRUE), sd(..., na.rm = TRUE)
    n <- bigstatsr::big_apply(corr_matrix,
        a.FUN = function(X, ind) {
            sum(!is.na(X[, ind, drop = FALSE])) # this works because TRUE = 1 and FALSE = 0 in R
        }, a.combine = "sum", block.size = 200
    )

    mean <- bigstatsr::big_apply(corr_matrix, a.FUN = function(X, ind) {
        sum(X[, ind, drop = FALSE], na.rm = TRUE)
    }, a.combine = "sum", block.size = 200) / n

    sd <- bigstatsr::big_apply(corr_matrix, a.FUN = function(X, ind) {
        sum((X[, ind, drop = FALSE] - mean)^2, na.rm = TRUE) / (n - 1)
    }, a.combine = "sum", block.size = 200) %>% sqrt()

    dstat <- list(
        mean = mean,
        sd = sd,
        sn = sd / mean # signal-to-noise ratio
    )
    return(dstat)
}

# --------------------------------------------------------------------------------------------- #
#' Calculates the mean matrix of a list of FBMs.
#' Every element is summed together and divided by the number of FBMs.
#'
#' This function takes a list of FBMs, and calculates the mean of every element in the FBM.
#' If the list is empty or the FBMs have different dimensions, it throws an error.
#'
#' @param fbmList A list of FBM objects to be added.
#' @return A new FBM object containing the sum of all FBMs in the list.
#' @examples
#' combined_matrix <- big_add_mat_list(list(fbm1, fbm2, fbm3))
#' @importFrom bigstatsr FBM
big_mat_list_mean <- function(fbmList) {
    # TODO: add input validations
    if (length(fbmList) <= 1) {
        stop("input list is empty or only contains one FBM")
    }

    # Check if all matrices in the list have the same dimensions
    if (!all(sapply(fbmList, function(x) identical(dim(fbmList[[1]]), dim(x))))) {
        stop("All matrices in the list have the same dimensions.")
    }
    n_col <- ncol(fbmList[[1]])
    n_row <- nrow(fbmList[[1]])
    mat_mean_zscore <- bigstatsr::FBM(n_row, n_col)

    bigstatsr::big_apply(mat_mean_zscore, a.FUN = function(X, ind) {
        X.sub <- X[, ind, drop = FALSE]
        X.sub <- Reduce(function(x, y) x + y[, ind, drop = FALSE], fbmList, init = X.sub)
        X[, ind] <- X.sub / length(fbmList)
        NULL
    }, a.combine = "c", block.size = 200)

    return(mat_mean_zscore)
}

# --------------------------------------------------------------------------------------------- #
#' Calculate statistical measures from a list of FBMs
#'
#' This function computes the mean, standard deviations, and coefficient of variation for each element
#' across a list of FBMs. The function first checks that the input list contains more than one FBM
#' and that all FBMs have the same dimensions. It calculates the z-score mean, standard deviation,
#' and coefficient of variation for each element across the FBMs in the list.
#'
#' @param fbmList A list of FBM objects from the bigstatsr package.
#'
#' @return A list containing three FBMs: mean_zscore, sds_zscore, and cv_zscore, which represent
#' the mean, standard deviation, and coefficient of variation z-scores of the elements across
#' the provided FBMs respectively.
#'
#' @importFrom bigstatsr FBM
#' @importFrom bigstatsr big_apply
#'
#' @examples
#' fbm_list <- list(fbm1, fbm2, fbm3)
#' stats_results <- get_list_stats(fbm_list)
#'
#' @details The function first calculates the mean for each element by summing up elements across
#' all FBMs and dividing by the number of FBMs. It then calculates the variance for each element,
#' takes the square root to get the standard deviation, and finally computes the coefficient of variation
#' for each element as the ratio of the standard deviation to the mean. The calculations avoid dividing by
#' zero by assigning NA to diagonal elements when computing the coefficient of variation, assuming
#' these might represent self-correlations.
#'
#' @seealso \code{\link[bigstatsr]{FBM}}, \code{\link[bigstatsr]{big_apply}}
#'
#' @export
get_list_stats <- function(fbmList) {
    if (length(fbmList) <= 1) {
        stop("input list is empty or only contains one FBM")
    }
    # Check if all matrices in the list have the same dimensions
    if (!all(sapply(fbmList, function(x) identical(dim(fbmList[[1]]), dim(x))))) {
        stop("All matrices in the list have the same dimensions.")
    }
    message("Calculating mean...")
    mat_mean_zscore <- big_mat_list_mean(fbmList)

    n_col <- ncol(fbmList[[1]])
    n_row <- nrow(fbmList[[1]])
    mat_sds_zscore <- bigstatsr::FBM(n_row, n_col)


    message("Calculating sds...")
    bigstatsr::big_apply(mat_sds_zscore, a.FUN = function(X, ind) {
        wrap_sds <- function(x, y) {
            x <- x + (y[, ind, drop = FALSE] - mat_mean_zscore[, ind, drop = FALSE])^2
        }
        X.sub <- X[, ind, drop = FALSE]
        X.sub <- Reduce(wrap_sds, fbmList, init = X.sub)
        X[, ind] <- sqrt(X.sub / (length(fbmList) - 1))
        NULL
    }, a.combine = "c", block.size = 200)

    mat_sn_zscore <- bigstatsr::FBM(n_row, n_col)
    diag(mat_sds_zscore) <- NA # to avoid dividing by zero ==> diagonal is 0 cause of perfect correlation ofc

    bigstatsr::big_apply(mat_sn_zscore, a.FUN = function(X, ind) {
        X[, ind] <- mat_mean_zscore[, ind, drop = FALSE] / mat_sds_zscore[, ind, drop = FALSE]
        NULL
    }, a.combine = "c", block.size = 200)
    res <- list(
        mean_zscore = mat_mean_zscore[],
        sds_zscore = mat_sds_zscore[],
        sn_zscore = mat_sn_zscore[] 
    )

    return(res)
}


# --------------------------------------------------------------------------------------------- #
select_genes <- function(
    lout,
    zscore_mean_threshold = 2,
    zscore_sn_threshold = 2,
    intersect_genes) {
    gene_ind <- which(
        abs(lout$list_stats$mean_zscore) > zscore_mean_threshold &
            abs(lout$list_stats$sn_zscore) > zscore_sn_threshold,
        arr.ind = TRUE
    )

    lout$gene_names <- intersect_genes[sort(unique(as.vector(gene_ind)))]
    return(lout)
}
