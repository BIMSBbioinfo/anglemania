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
#' @param extrema double. Fraction of the correlation
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
                          extrema) {
    if (is.null(rownames(x_mat))) {
        stop("Input matrix lacks rownames. Stopping")
    }
    if (length(extrema) > 2 || length(extrema) < 1) {
        stop("extrema should be numeric of length either 1 or 2")
    }
    gene_names <- rownames(x_mat)
    n_seacells <- ncol(x_mat)
    # compute correlation matrix
    
    x_mat_ang <- big_extract_corr(
        x_mat %>% as.matrix()
    )

    # Estimate critical angles ==> Using a Gaussian fit.

    l_angles <- big_estimate_critical_angles(x_mat_ang,
        s_dims = n_seacells,
        extrema = extrema
    )

    # TODO: make something like a progress bar

    results_dt <- big_apply(x_mat_ang,
        a.FUN = function(x_mat_ang, ind) {

            # create subset matrix
            X.sub <- x_mat_ang[, ind, drop = FALSE]


            cols_matrix <- matrix(rep(ind, each = nrow(x_mat_ang)), nrow = nrow(x_mat_ang), ncol = length(ind), byrow = FALSE)
            # this is basically equivalent to cols(matrix) but the problem with the big_apply method is that
            # for each block the function is applied to, the columns are  1:500 for every block. (if blocksize = 500)
            # row(X.sub) <= cols_matrix is equivalent to upper.tri(X.sub)
            sharp <- which(X.sub <= min(l_angles$critical_angles) & (row(X.sub) <= cols_matrix), arr.ind = TRUE) # get indices of the sharp angles
            blunt <- which(X.sub >= max(l_angles$critical_angles) & (row(X.sub) <= cols_matrix), arr.ind = TRUE) # get indices of the blunt angles

            sharp_values <- X.sub[sharp]
            blunt_values <- X.sub[blunt]

            # Adjust indices in the results to reflect the original matrix
            # and convert rownumber/colnumber to gene names
            sharp[, 2] <- gene_names[ind[sharp[, 2]]]
            sharp[, 1] <- gene_names[sharp[, 1] %>% as.integer()] # somehow the indices are not integers, but only in the second column? lel

            blunt[, 2] <- gene_names[ind[blunt[, 2]]]
            blunt[, 1] <- gene_names[blunt[, 1] %>% as.integer()]


            # Combine indices and values
            sharp_data <- data.table(
                gene_a = sharp[, 1],
                gene_b = sharp[, 2],
                angle = sharp_values,
                angle_type = "sharp",
                signf = 1
            )

            blunt_data <- data.table(
                gene_a = blunt[, 1],
                gene_b = blunt[, 2],
                angle = blunt_values,
                angle_type = "blunt",
                signf = 1
            )

            # Combine sharp and blunt data
            rbind(sharp_data, blunt_data)
        }, a.combine = "rbind", block.size = 1000
    )
    # create a column in the data.table specifying the name of the dataset
    # COMMENT: the name of the dataset is currently specified as the name, the mclapply function iterates over
    # FIXME: where dataset specification ==> name(list_element)
    results_dt$dataset <- name

    results_list <- list(results = results_dt, l_angles = l_angles)
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
    return(results_list)
}


