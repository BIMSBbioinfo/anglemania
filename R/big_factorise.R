#' Factorise angle matrices
#'
#' @description
#' Generates two a one-hot encoded square matrixces
#' recording the presence/absence (1/0) of sharp and blunt
#' critical angles.
#'
#' @details
#' *factorise* extracts angles between genes, estimates
#' critical angles by approximation of the angle distribution,
#' and records angles passing the critical threshold into
#' a sparse matrix.
#'
#' @param x_mat Matrix. Contains normalised and scaled gene
#'   expression.
#' @param name character. The name of the dataset.
#' @param fdr_threshold double. Fraction of the correlation
#'   to be cut from both sides of an approximated angle
#'   distribution.
#' @return list. First two elements are sparse matrices
#'   containing the significant sharp and blunt angles
#'   between genes. Third and fourth elements are lists
#'   with angles statistics and paths to the values of
#'   critical mangles.
#' @export big_factorise
big_factorise <- function(x_mat, # nolint
                          name,
                          fdr_threshold = 0.005) {
    if (is.null(rownames(x_mat))) {
        stop("Input matrix lacks rownames. Stopping")
    }
    # if (length(extrema) > 2 || length(extrema) < 1) {
    #     stop("extrema should be numeric of length either 1 or 2")
    # }
    gene_names <- rownames(x_mat)
    n_seacells <- ncol(x_mat)
    # compute correlation matrix
    
    x_mat_corr <- big_extract_corr(
        x_mat %>% as.matrix()
    )


    # TODO: make something like a progress bar

    results_dt <- big_apply(x_mat_corr,
        a.FUN = function(x_mat_corr, ind) {

            # create subset matrix
            X.sub <- x_mat_corr[, ind, drop = FALSE]

            cols_matrix <- matrix(rep(ind, each = nrow(x_mat_corr)), nrow = nrow(x_mat_corr), ncol = length(ind), byrow = FALSE)
            # this is basically equivalent to cols(matrix) but the problem with the big_apply method is that
            # for each block the function is applied to, the columns are  1:500 for every block. (if blocksize = 500)
            # row(X.sub) <= cols_matrix is equivalent to upper.tri(X.sub)
            # select the top genes with the lowest q-value
            keep <- which(X.sub <= fdr_threshold & (row(X.sub) <= cols_matrix), arr.ind = TRUE) # get indices of the sharp angles

            keep_values <- X.sub[keep]
            # print(keep_values[1:5])
            # Adjust indices in the results to reflect the original matrix
            # and convert rownumber/colnumber to gene names
            keep[, 2] <- gene_names[ind[keep[, 2]]]
            # print(keep[, 2][1:5])
            keep[, 1] <- gene_names[keep[, 1] %>% as.integer()] # somehow the indices are not integers, but only in the second column? lel

            # print(keep[, 1][1:5])
            # Combine indices and values
            keep_data <- data.table(
                gene_a = keep[, 1],
                gene_b = keep[, 2],
                q_value = keep_values
            )
        }, a.combine = "rbind", block.size = 1000
    )
    # create a column in the data.table specifying the name of the dataset
    # COMMENT: the name of the dataset is currently specified as the name, the mclapply function iterates over
    # FIXME: where dataset specification ==> name(list_element)
    results_dt$dataset <- name
    
    # COMMENT: could now either work with the absolute angles or binarize like before
    # The Results list includes the following elements:
    # 1. results_dt ==> a data.table with the gene pairs with sharp and blunt angles and their corresponding angle
    # Sample data.table representation:
    # ---------------------------------------------
    # |  gene_a  |  gene_b  |  angle   | angle_type |  dataset  |
    # ---------------------------------------------
    # |  gene1   |  gene2   |  12.33   |   sharp    | dataset_X |
    # |  gene3   |  gene4   |  8.74    |   sharp    | dataset_X |
    # |  gene5   |  gene6   |  15.64   |   sharp    | dataset_X |
    # ---------------------------------------------
    # 2. l_angles ==> a list with information about the calculation of the critical angles
    # used for the filtering of the sharp and blunt angles

    invisible(gc())
    return(results_dt)
}


