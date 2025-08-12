sce_raw <- sce_example()

test_that("sparse_to_fbm converts sparse matrix to FBM correctly", {
    # Generate a random sparse matrix using the Matrix package
    s_mat <- Matrix::rsparsematrix(nrow = 10, ncol = 5, density = 0.3)

    # Convert the sparse matrix to an FBM using your function
    fbm_mat <- sparse_to_fbm(s_mat)

    # Check that the result is an FBM object
    expect_s4_class(fbm_mat, "FBM")

    # Verify that the dimensions match
    expect_equal(dim(fbm_mat), dim(s_mat))

    # Compare the content of the FBM with the original sparse matrix
    # Coerce the FBM to a regular matrix for comparison
    fbm_data <- fbm_mat[]
    s_mat_dense <- as.matrix(s_mat)

    expect_equal(fbm_data, s_mat_dense)
})

# File: tests/testthat/test-get_dstat.R

test_that("get_dstat computes statistics correctly for an FBM", {
    # Create a small test matrix with known values
    # Important to know that we deal with symmetric matrices!!!
    mat <- matrix(
        c(
            NA,
            0.1,
            -0.2,
            -1,
            0.1,
            -0.5,
            0.7,
            -0.3,
            0.2,
            -0.1,
            0.1,
            NA,
            -1,
            0.1,
            -0.3,
            -0.2,
            0.5,
            -0.1,
            -0.2,
            0.3,
            -0.2,
            -1,
            NA,
            0.2,
            -0.1,
            0.4,
            -0.3,
            -0.2,
            0.1,
            -0.5,
            -1,
            0.1,
            0.2,
            NA,
            0.3,
            -0.5,
            0.1,
            -0.2,
            0.7,
            0.3,
            0.1,
            -0.3,
            -0.1,
            0.3,
            NA,
            -0.2,
            -0.5,
            0.1,
            -0.7,
            0.2,
            -0.5,
            -0.2,
            0.4,
            -0.5,
            -0.2,
            NA,
            0.3,
            0.2,
            -0.1,
            0.4,
            0.7,
            0.5,
            -0.3,
            0.1,
            -0.5,
            0.3,
            NA,
            -0.2,
            -0.1,
            -0.3,
            -0.3,
            -0.1,
            -0.2,
            -0.2,
            0.1,
            0.2,
            -0.2,
            NA,
            0.5,
            -0.4,
            0.2,
            -0.2,
            0.1,
            0.7,
            -0.7,
            -0.1,
            -0.1,
            0.5,
            NA,
            -0.3,
            -0.1,
            0.3,
            -0.5,
            0.3,
            0.2,
            0.4,
            -0.3,
            -0.4,
            -0.3,
            NA
        ),
        nrow = 10,
        ncol = 10,
        byrow = TRUE
    )

    # Convert the matrix to an FBM
    fbm_mat <- bigstatsr::FBM(
        nrow = nrow(mat),
        ncol = ncol(mat),
        init = mat
    )

    # Compute the statistics using get_dstat
    result <- get_dstat(fbm_mat)

    # Manually compute the expected statistics
    n_expected <- length(mat[!is.na(mat)])
    # this n is used in anglemania because
    # we use the entire matrix not only the upper.tri
    mean_expected <- Matrix::colMeans(mat, na.rm = TRUE)
    var_expected <- matrixStats::colVars(mat, na.rm = TRUE)
    sd_expected <- matrixStats::colSds(mat, na.rm = TRUE)
    min_expected <- matrixStats::colMins(mat, na.rm = TRUE)
    max_expected <- matrixStats::colMaxs(mat, na.rm = TRUE)

    # Check that the computed values match the expected values
    expect_equal(result$mean, mean_expected)
    expect_equal(result$var, var_expected)
    expect_equal(result$sd, sd_expected)
    expect_equal(result$min, min_expected)
    expect_equal(result$max, max_expected)
})

test_that("get_dstat throws an error when corr_matrix is not an FBM", {
    # Create a regular matrix
    regular_matrix <- matrix(1:9, nrow = 3)

    # Expect an error when passing a non-FBM object
    expect_error(
        get_dstat(regular_matrix),
        "corr_matrix has to be a FBM object"
    )
})

test_that("big_mat_list_mean computes the weighted mean correctly", {
    library(S4Vectors)
    sce <- sce_raw
    batch_key <- "batch"
    sce <- anglemania(sce, batch_key = batch_key)
    # Compute the weighted mean using big_mat_list_mean
    weights <- setNames(
        S4Vectors::metadata(sce)$anglemania$weights$weight,
        S4Vectors::metadata(sce)$anglemania$weights$anglemania_batch
    )
    result_fbm <- big_mat_list_mean(
        metadata(sce)$anglemania$matrix_list,
        weights = weights
    )

    # Extract the result as a regular matrix ==> only use a subset
    result_matrix <- result_fbm[1:100, 1:100]

    # Manually compute the expected weighted mean ==> only use a subset
    expected_matrix <- metadata(sce)$anglemania$matrix_list$batch1[
        1:100,
        1:100
    ] *
        0.5 +
        metadata(sce)$anglemania$matrix_list$batch2[1:100, 1:100] *
            0.5

    # Compare the result with the expected matrix
    expect_equal(result_matrix, expected_matrix)
})


test_that("big_mat_list_mean throws an error
when matrices have different dimensions", {
    # Create FBMs with different dimensions
    mat1 <- matrix(1:9, nrow = 3)
    mat2 <- matrix(1:8, nrow = 4) # Different dimensions

    fbm1 <- bigstatsr::FBM(
        nrow = nrow(mat1),
        ncol = ncol(mat1),
        init = mat1
    )
    fbm2 <- bigstatsr::FBM(
        nrow = nrow(mat2),
        ncol = ncol(mat2),
        init = mat2
    )

    # Create weights
    weights <- c(batch1 = 0.5, batch2 = 0.5)

    # Create the list of FBMs
    fbm_list <- list(batch1 = fbm1, batch2 = fbm2)

    # Expect an error
    expect_error(
        big_mat_list_mean(fbm_list, weights = weights),
        "All matrices in the list must have the same dimensions."
    )
})


test_that("select_genes selects genes correctly based on thresholds", {
    library(SingleCellExperiment)
    sce <- sce_raw
    # Create minimal matrices for mean_zscore and sn_zscore
    set.seed(123) # For reproducibility
    mean_zscore_matrix <- bigstatsr::FBM(
        nrow = 20,
        ncol = 20,
        init = rnorm(400, mean = 0, sd = 3)
    )
    sn_zscore_matrix <- bigstatsr::FBM(
        nrow = 20,
        ncol = 20,
        init = rnorm(400, mean = 0, sd = 1)
    )

    # Define intersect_genes
    gene_names <- paste0("gene", 1:20)

    metadata(sce)$anglemania$list_stats <- list(
        mean_zscore = mean_zscore_matrix,
        sn_zscore = sn_zscore_matrix
    )
    metadata(sce)$anglemania$intersect_genes <- gene_names
    metadata(sce)$anglemania$prefiltered_df <- prefilter_angl(
        snr_zscore_matrix = sn_zscore_matrix,
        mean_zscore_matrix = mean_zscore_matrix,
        zscore_mean_threshold = 1.0,
        zscore_sn_threshold = 1.0
    )

    # Call select_genes with thresholds
    sce <- select_genes(
        sce = sce,
        zscore_mean_threshold = 2.0,
        zscore_sn_threshold = 2.0,
        max_n_genes = 3
    )

    # Check that the integration_genes slot is updated correctly
    selected_genes <- get_anglemania_genes(sce)
    selected_info <- get_anglemania_stats_df(sce)

    # Expected gene indices (since we have 3 genes, and thresholds are 2.0)
    # We need to find gene pairs where both mean_zscore
    #  and sn_zscore exceed thresholds
    # For upper triangular matrix positions

    # Calculate expected indices from the prefiltered_df
    expected_indices <- which(
        upper.tri(sn_zscore_matrix[]) &
            (sn_zscore_matrix[] >= 2.0) &
            (abs(mean_zscore_matrix[]) >= 2.0),
        arr.ind = TRUE
    )

    # Extract expected gene pairs
    expected_info_df <- data.frame(
        geneA = gene_names[expected_indices[, 1]],
        geneB = gene_names[expected_indices[, 2]],
        sn_zscore = sn_zscore_matrix[expected_indices],
        mean_zscore = mean_zscore_matrix[expected_indices]
    )

    # Order by absolute zscore
    expected_info_df <- expected_info_df[
        order(abs(expected_info_df$mean_zscore), decreasing = TRUE),
    ]
    rownames(expected_info_df) <- NULL

    # Limit to max_n_genes
    expected_selected_genes <- unique(as.vector(rbind(
        expected_info_df$geneA,
        expected_info_df$geneB
    )))[1:3]

    # # Expected selected genes
    # expected_selected_genes <- unique(c(
    #   expected_gene_pairs
    # ))

    # Compare the selected genes
    testthat::expect_equal(selected_genes, expected_selected_genes)

    # Compare the selected info
    testthat::expect_equal(selected_info, expected_info_df)
})
