# --------------------------------------------------------------------------------------------- #
#' Permute counts in a Seurat object
#'
#' This function permutes the counts in the assay data of a Seurat object either by columns or rows.
#' @param seu A Seurat object.
#' @param which A character string specifying whether to permute "columns" or "rows".
#' @return A new Seurat object with permuted count data.
#' @examples
#' seu_permuted = permute_counts(seu, which = "columns")
permute_counts = function(seu, which = "columns", seed=1) {
    
    cnts = GetAssayData(seu, "RNA", "counts")
    if(which == "columns"){
        cnts = apply(cnts, 2, sample)
    } else if(which == "rows"){
        cnts = t(apply(cnts, 1, sample))
    }
    rownames(cnts) = rownames(seu)
    CreateSeuratObject(
        counts = cnts, 
        meta.data = seu@meta.data
    )
}

# --------------------------------------------------------------------------------------------- #
seurat_list_permute_counts = function(seurat_list, which = "columns", seed = 1){
    set.seed(seed)
    lapply(seurat_list, permute_counts, which=which)
}

# --------------------------------------------------------------------------------------------- #
#' Remove lower triangle of a matrix
#'
#' This function sets the lower triangular part of a matrix, including the diagonal, to NA.
#' @param mat A matrix.
#' @return A matrix with its lower triangular part set to NA.
#' @examples
#' mat = matrix(1:9, 3, 3)
#' mat_no_lower = remove_lower_tri_mat(mat)
remove_lower_tri_mat = function(mat) {
    mat[lower.tri(mat, diag = TRUE)] = NA
    return(mat)
}


# --------------------------------------------------------------------------------------------- #
#' Calculate mean and standard deviation of the upper triangle of a matrix
#'
#' This function computes the mean and standard deviation of the upper triangular part of the matrix, excluding NAs.
#' @param mat A matrix.
#' @return A data frame with the mean and standard deviation.
#' @examples
#' mat = matrix(1:9, 3, 3)
#' stats = matrix_mean_sd(mat)
matrix_mean_sd = function(mat) {
    upper_tri_mat = remove_lower_tri_mat(mat)
    return(data.frame(
        mean = mean(upper_tri_mat, na.rm = TRUE),
        sds  = sd(upper_tri_mat, na.rm = TRUE)
    ))
}

# --------------------------------------------------------------------------------------------- #
#' Calculate mean and standard deviation for a list of matrices
#'
#' This function applies a matrix mean and standard deviation calculation across a list of matrices and combines the results.
#' @param matrix_list A list of matrices.
#' @return A data frame of combined results from each matrix.
#' @examples
#' matrix_list_stats = matrix_list_mean_sds(list(matrix(1:4, 2, 2), matrix(5:8, 2, 2)))
matrix_list_mean_sds = function(matrix_list) {
    lstat = lapply(matrix_list, function(x) {
        matrix_mean_sd(x)
    })
    dstat = do.call(rbind, lstat)
    dstat
}

# --------------------------------------------------------------------------------------------- #
#' Standardize matrix using z-score
#'
#' This function standardizes the elements of a matrix based on provided mean and standard deviation.
#' @param mat A matrix to be standardized.
#' @param dstat A data frame with 'mean' and 'sd' columns used for standardization.
#' @return The standardized matrix.
#' @examples
#' standardized_mat = matrix_z_score(mat, data.frame(mean = 5, sd = 2))
matrix_z_score = function(mat, dstat) {
   (mat - dstat$mean) / dstat$sd
}

# --------------------------------------------------------------------------------------------- #
#' Compute z-scores for a list of matrices
#'
#' This function applies z-score standardization across a list of matrices based on provided statistical data for each matrix.
#' @param matrix_list A list of matrices.
#' @param dstat A data frame or list containing the mean and standard deviation for each matrix, used for standardization.
#' @return A list of matrices standardized by z-score.
#' @examples
#' zscores = matrix_list_z_score(list_of_matrices, stats_data_frame)
matrix_list_z_score = function(matrix_list, dstat) {
    matrix_list_zscore = list()
    for(i in seq(matrix_list)) {
        matrix_list_zscore[[as.character(i)]] = matrix_z_score(matrix_list[[i]], dstat[i,])
    }
    return(matrix_list_zscore)
}

# --------------------------------------------------------------------------------------------- #
#' Calculate mean, standard deviation, and coefficient of variation of z-scores for a list of matrices
#'
#' This function calculates the mean, standard deviation, and coefficient of variation of z-scores across a list of matrices after setting NAs to zero in the lower triangular part.
#' @param matrix_list_zscore A list of matrices with z-scores.
#' @return A list containing the mean z-score, standard deviation of z-scores, and coefficient of variation of z-scores.
#' @examples
#' summary_stats = matrix_list_z_score_mean_sd(list_of_zscore_matrices)
matrix_list_z_score_mean_sd = function(matrix_list_zscore) {

    message("remove lower tri ..")
    matrix_list_zscore = lapply(matrix_list_zscore, remove_lower_tri_mat)
    
    message("remove NA ..")
    matrix_list_zscore = lapply(matrix_list_zscore, function(mat) {
        mat[is.na(mat)] = 0
        mat
    })

    message("mean ..")
    mat_mean_zscore = matrix_list_zscore[[1]]
    for(i in 2:length(matrix_list_zscore)) {
        mat_mean_zscore = mat_mean_zscore + matrix_list_zscore[[i]]
    }
    mat_mean_zscore = mat_mean_zscore / length(matrix_list_zscore)

    message("sd ..")
    mat_sds_zscore = (matrix_list_zscore[[1]] - mat_mean_zscore)^2
    for(i in 2:length(matrix_list_zscore)) {
        mat_sds_zscore = mat_sds_zscore + (matrix_list_zscore[[i]] - mat_mean_zscore)^2
    }
    mat_sds_zscore = sqrt(mat_sds_zscore / (length(matrix_list_zscore) - 1))

    message("returning ..")
    return(list(
        mean_zscore = mat_mean_zscore,
        sds_zscore  = mat_sds_zscore,
        cv_zscore   = mat_mean_zscore / mat_sds_zscore
    ))
}

# --------------------------------------------------------------------------------------------- #
#' Calculate angles for scale data matrices from a list of Seurat objects
#'
#' This function extracts scale data matrices from each Seurat object in a list and then calculates angles using a specified number of threads.
#' @param seurat_list A list of Seurat objects.
#' @return A list of results from the angle calculation for each Seurat object's scale data matrix.
#' @examples
#' angle_results = calculate_angles_seurat_list(list_of_seurats)
calculate_angles_seurat_list = function(seurat_list) {
    lapply(seurat_list, function(seu) {
        if(!"scale.data" %in% Layers(seu))
            stop("scale.data not in seu")
        x_mat = GetAssayData(seu, "RNA", "scale.data")
        extract_angles(x_mat, n_threads = 16)
    })
}


# --------------------------------------------------------------------------------------------- #
select_genes = function(
    lout,
    zscore_mean_threshold = 2,
    zscore_cv_threshold   = 2

){
    gene_ind = which(
        abs(lout$l_zscore_mean_sd$mean_zscore)  > zscore_mean_threshold & 
        abs(lout$l_zscore_mean_sd$cv_zscore) > zscore_cv_threshold, 
        arr.ind=TRUE
    )
    lout$gene_names = rownames(lout$l_zscore_mean_sd$mean_zscore)[sort(unique(as.vector(gene_ind)))]
    return(lout)

}