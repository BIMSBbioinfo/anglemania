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
big_mat_list_mean <- function(anglem_object) {
    # TODO: add input validations 
    if (class(anglem_object) != "anglem") {
        stop("anglem_object needs to be a anglem object")
    }
    # Check if all matrices in the list have the same dimensions
    if (!all(sapply(matrix_list(anglem_object), function(x) identical(dim(matrix_list(anglem_object)[[1]]), dim(x))))) {
        stop("All matrices in the list have the same dimensions.")
    }

    n_col <- ncol(matrix_list(anglem_object)[[1]])
    n_row <- nrow(matrix_list(anglem_object)[[1]])
    mat_mean_zscore <- bigstatsr::FBM(n_row, n_col)

    bigstatsr::big_apply(mat_mean_zscore, a.FUN = function(X, ind) {
        X.sub <- X[, ind, drop = FALSE]
        wrap_mean <- function(final_mat, batch) {
            batch_mat <- matrix_list(anglem_object)[[batch]][, ind, drop = FALSE] * anglem_object@weights[batch]
            final_mat <- final_mat + batch_mat
        }
        #COMMENT: So currently I run this function in the reduce statement on the names of the weights vector
        X.sub <- Reduce(wrap_mean, names(anglem_object@weights), init = X.sub)
        X[, ind] <- X.sub # / length(matrix_list(anglem_object))
        #COMMENT: 08/07/2024 now that the individual zscore matrices have been scaled by the weight, which is adjusted so that sum(weight)==1 we do not have to
        #COMMENT: divide by the number of matrices 
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
get_list_stats <- function(anglem_object) {
    # TODO: change for anglem object    
    if (class(anglem_object) != "anglem") {
        stop("anglem_object needs to be a anglem object")
    }
    # Check if all matrices in the list have the same dimensions
    if (!all(sapply(matrix_list(anglem_object), function(x) identical(dim(matrix_list(anglem_object)[[1]]), dim(x))))) {
        stop("All matrices in the list need to have the same dimensions.")
    }
    message("Weighting matrix_list...")
    # Weight all matrices
    # tmp_names <- names(anglem_object@weights)
    # matrix_list(anglem_object) <- pbapply::pblapply(tmp_names, function(batch_name) {
    #     print(batch_name)
    #     bigstatsr::big_apply(matrix_list(anglem_object)[[batch_name]], function(X, ind) {
    #         X[, ind] <- X[, ind] * anglem_object@weights[batch_name]
    #         NULL
    #     }, a.combine = "c", block.size = 200)
    #     return(matrix_list(anglem_object)[[batch_name]])
    # })
    # names(matrix_list) <- tmp_names
    message("Calculating mean...")
    mat_mean_zscore <- big_mat_list_mean(anglem_object)
    
    n_col <- ncol(matrix_list(anglem_object)[[1]])
    n_row <- nrow(matrix_list(anglem_object)[[1]])
    mat_sds_zscore <- bigstatsr::FBM(n_row, n_col)


    message("Calculating sds...")
    # TODO: add weighting
    bigstatsr::big_apply(mat_sds_zscore, a.FUN = function(X, ind) {
        wrap_sds <- function(final_mat, batch) {
            batch_mat <- matrix_list(anglem_object)[[batch]]
            final_mat <- final_mat + (batch_mat[, ind, drop = FALSE] - mat_mean_zscore[, ind, drop = FALSE])^2 * anglem_object@weights[batch]
        }
        X.sub <- X[, ind, drop = FALSE]
        X.sub <- Reduce(wrap_sds, names(anglem_object@weights), init = X.sub)
        X[, ind] <- sqrt(X.sub) # / length(matrix_list(anglem_object))
        #COMMENT: 08/07/2024 now that the individual zscore matrices have been scaled by the weight, which is adjusted so that sum(weight)==1 we do not have to
        #COMMENT: divide by the number of matrices
        NULL
    }, a.combine = "c", block.size = 200)

    mat_sn_zscore <- bigstatsr::FBM(n_row, n_col)
    diag(mat_sds_zscore) <- NA # to avoid dividing by zero ==> diagonal is 0 cause of perfect correlation ofc

    bigstatsr::big_apply(mat_sn_zscore, a.FUN = function(X, ind) {
        X[, ind] <- abs(mat_mean_zscore[, ind, drop = FALSE]) / mat_sds_zscore[, ind, drop = FALSE]
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
extract_rows_for_unique_genes <- function(dt, max_n_genes) {
    unique_genes <- numeric()
    for (i in 1:nrow(dt)) {
        # Add new unique genes from both columns of the current row
        current_genes <- c(dt$geneA[i], dt$geneB[i])
        unique_genes <- unique(c(unique_genes, current_genes))

        # Stop if we have accumulated max_n_genes or more unique genes
        if (length(unique_genes) >= max_n_genes) {
            break
        }
    }
    unique_genes <- sort(unique_genes)
    return(unique_genes)
}

select_genes <- function(
    anglem_object,
    zscore_mean_threshold = 2,
    zscore_sn_threshold = 2,
    max_n_genes = 2000) {

    if (class(anglem_object) != "anglem") {
        stop("anglem_object needs to be a anglem object")
    }

    gene_ind <- which(
        upper.tri(list_stats(anglem_object)$sn_zscore) & 
        (list_stats(anglem_object)$sn_zscore >= zscore_sn_threshold) & # signal-to-noise ratio is above 0 anyways, no need to use absolute
        (abs(list_stats(anglem_object)$mean_zscore) >= zscore_mean_threshold), # for mean zscore, use absolute value
        arr.ind = TRUE
    )

    #COMMENT: 08/07/2024 current approach:
    # If no genes passed the cutoff, increase the threshold
    if (nrow(gene_ind) == 0) {
        message("No genes passed the cutoff.")
        quantile90mean <- quantile(abs(list_stats(anglem_object)$mean_zscore), 0.95, na.rm = TRUE)
        quantile90sn <- quantile(list_stats(anglem_object)$sn_zscore, 0.95, na.rm = TRUE)
        if (quantile90mean < zscore_mean_threshold) {
            zscore_mean_threshold <- quantile90mean
            message(paste0("zscore_mean_threshold is lower than the 95% quantile of the absolute mean zscores. Setting the zscore_mean_threshold to:\t", zscore_mean_threshold))
        }
        if (quantile90sn < zscore_sn_threshold) {
            zscore_sn_threshold <- quantile90sn
            message(paste0("zscore_sn_threshold is higher than 95% quantile. Setting the zscore_sn_threshold to:\t", round(zscore_sn_threshold, 2)))
        }
        print(paste0("zscore_mean_threshold:", zscore_mean_threshold, " zscore_sn_threshold:", round(zscore_sn_threshold, 2)))

        # Re-run the selection
        gene_ind <- which(
            upper.tri(list_stats(anglem_object)$sn_zscore) & 
            (list_stats(anglem_object)$sn_zscore >= zscore_sn_threshold) & # signal-to-noise ratio is above 0 anyways, no need to use absolute
            (abs(list_stats(anglem_object)$mean_zscore) >= zscore_mean_threshold), # for mean zscore, use absolute value
            arr.ind = TRUE
        )

        message("If desired, you can re-run the selection of genes with a lower zscore_mean_threshold and/or zscore_sn_threshold by using the 'select_genes' function. e.g.: anglem_object <- select_genes(anglem_object, zscore_mean_threshold = 1, zscore_sn_threshold = 1, max_n_genes = 2000)")
        message("Please inspect integration_genes(anglem_object)$info for info on the scores of the selected gene pairs.")
    }

    top_n <- data.frame(
        geneA = gene_ind[, 1],
        geneB = gene_ind[, 2],
        zscore = list_stats(anglem_object)$mean_zscore[gene_ind],
        snscore = list_stats(anglem_object)$sn_zscore[gene_ind]
    )


    #order df
    top_n <- top_n[order(top_n$zscore, decreasing = TRUE), ]
    anglem_object@integration_genes$info <- top_n
    top_n <- extract_rows_for_unique_genes(top_n, max_n_genes)
    anglem_object@integration_genes$genes <- intersect_genes(anglem_object)[top_n]
    return(anglem_object)
}
