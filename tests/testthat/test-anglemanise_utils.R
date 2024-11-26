# load example seurat object
# pbmc_scmall does not have any "batches" from multiple experiments
#   but we'll just treat it like the "groups" column are the batches
se <- SeuratObject::pbmc_small
se_raw <- se

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
      NA, 0.1, -0.2, -1, 
      0.1, NA, -1, -0.1, 
      -0.2, -1, NA, 0.2, 
      -1, -0.1, 0.2, NA
    ),
    nrow = 4,
    ncol = 4,
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
  mean_expected <- mean(mat[upper.tri(mat, diag = FALSE)], na.rm = TRUE)
  var_expected <- var(mat[upper.tri(mat, diag = FALSE)], na.rm = TRUE)
  sd_expected <- sd(mat[upper.tri(mat, diag = FALSE)], na.rm = TRUE)
  min_expected <- min(mat[upper.tri(mat, diag = FALSE)], na.rm = TRUE)
  max_expected <- max(mat[upper.tri(mat, diag = FALSE)], na.rm = TRUE)
  sn_expected <- sd_expected / mean_expected

  # Check that the computed values match the expected values
  expect_equal(result$mean, mean_expected)
  expect_equal(result$var, var_expected)
  expect_equal(result$sd, sd_expected)
  expect_equal(result$sn, sn_expected)
  expect_equal(result$min, min_expected)
  expect_equal(result$max, max_expected)
})

test_that("get_dstat throws an error when corr_matrix is not an FBM", {
  # Create a regular matrix
  regular_matrix <- matrix(1:9, nrow = 3)

  # Expect an error when passing a non-FBM object
  expect_error(get_dstat(regular_matrix), "corr_matrix has to be a FBM object")
})

test_that("big_mat_list_mean computes the weighted mean correctly", {
  se <- se_raw
  batch_key <- "groups"
  angl <- create_anglemania_object(
    se,
    batch_key = batch_key
  )
  angl <- anglemania(angl)
  # Compute the weighted mean using big_mat_list_mean
  result_fbm <- big_mat_list_mean(angl)

  # Extract the result as a regular matrix ==> only use a subset
  result_matrix <- result_fbm[1:100, 1:100]

  # Manually compute the expected weighted mean ==> only use a subset
  expected_matrix <- angl@matrix_list$g1[1:100, 1:100] * 0.5 +
    angl@matrix_list$g2[1:100, 1:100] * 0.5


  # Compare the result with the expected matrix
  expect_equal(result_matrix, expected_matrix)
})


test_that("big_mat_list_mean throws an error
when matrices have different dimensions", {
  # Create FBMs with different dimensions
  mat1 <- matrix(1:9, nrow = 3)
  mat2 <- matrix(1:8, nrow = 4) # Different dimensions

  fbm1 <- bigstatsr::FBM(nrow = nrow(mat1), ncol = ncol(mat1), init = mat1)
  fbm2 <- bigstatsr::FBM(nrow = nrow(mat2), ncol = ncol(mat2), init = mat2)

  # Create weights
  weights <- c(batch1 = 0.5, batch2 = 0.5)

  # Create the list of FBMs
  fbm_list <- list(batch1 = fbm1, batch2 = fbm2)

  # Construct the anglemania_object
  angl <- new(
    "anglemania_object",
    weights = weights,
    matrix_list = fbm_list
  )

  # Expect an error
  expect_error(
    big_mat_list_mean(angl),
    "All matrices in the list must have the same dimensions."
  )
})


test_that("select_genes selects genes correctly based on thresholds", {
  # Create minimal matrices for mean_zscore and sn_zscore
  set.seed(123) # For reproducibility
  mean_zscore_matrix <- matrix(
    rnorm(400, mean = 0, sd = 3), # Mean 0, SD 3
    nrow = 20,
    ncol = 20
  )

  sn_zscore_matrix <- matrix(
    rnorm(400, mean = 2, sd = 1), # Mean 2, SD 1
    nrow = 20,
    ncol = 20
  )

  # Define intersect_genes
  gene_names <- paste0("gene", 1:20)

  # Create an anglemania_object with the required slots
  test_angl <- new(
    "anglemania_object",
    list_stats = list(
      mean_zscore = mean_zscore_matrix,
      sn_zscore = sn_zscore_matrix
    ),
    intersect_genes = gene_names,
    integration_genes = list(
      info = data.frame(),
      genes = character()
    )
  )

  # Call select_genes with thresholds
  updated_object <- select_genes(
    angl = test_angl,
    zscore_mean_threshold = 2.0,
    zscore_sn_threshold = 2.0,
    max_n_genes = 3
  )

  # Check that the integration_genes slot is updated correctly
  selected_genes <- updated_object@integration_genes$genes
  selected_info <- updated_object@integration_genes$info

  # Expected gene indices (since we have 3 genes, and thresholds are 2.0)
  # We need to find gene pairs where both mean_zscore 
  #  and sn_zscore exceed thresholds
  # For upper triangular matrix positions

  # Calculate expected indices
  expected_indices <- which(
    upper.tri(sn_zscore_matrix) &
      (sn_zscore_matrix >= 2.0) &
      (abs(mean_zscore_matrix) >= 2.0),
    arr.ind = TRUE
  )

  # Extract expected gene pairs
  expected_info_df <- data.frame(
    geneA = expected_indices[, 1],
    geneB = expected_indices[, 2],
    zscore = mean_zscore_matrix[expected_indices],
    snscore = sn_zscore_matrix[expected_indices]
  )

  # Order by absolute zscore
  expected_info_df <- expected_info_df[
    order(abs(expected_info_df$zscore), decreasing = TRUE),
  ]

  # Limit to max_n_genes
  expected_gene_pairs <- unique(as.vector(rbind(
    expected_info_df[seq_len(2), ]$geneA,
    expected_info_df[seq_len(2), ]$geneB
  )))[1:3]

  # Expected selected genes
  expected_selected_genes <- unique(c(
    gene_names[expected_gene_pairs]
  ))

  # Compare the selected genes
  expect_equal(selected_genes, expected_selected_genes)

  # Compare the info data frame
  expect_equal(selected_info, expected_info_df)
})


test_that("select_genes adjusts thresholds when no genes pass the cutoff with larger matrices", {
  # Use higher thresholds so that no genes pass
  # Create minimal matrices for mean_zscore and sn_zscore
  set.seed(123) # For reproducibility
  mean_zscore_matrix <- matrix(
    rnorm(1600, mean = 0, sd = 3), # Mean 0, SD 3
    nrow = 40,
    ncol = 40
  )

  sn_zscore_matrix <- matrix(
    rnorm(1600, mean = 2, sd = 1), # Mean 2, SD 1
    nrow = 40,
    ncol = 40
  )

  # Define intersect_genes
  gene_names <- paste0("gene", 1:40)

  # Create an anglemania_object with the required slots
  test_angl <- new(
    "anglemania_object",
    list_stats = list(
      mean_zscore = mean_zscore_matrix,
      sn_zscore = sn_zscore_matrix
    ),
    intersect_genes = gene_names,
    integration_genes = list(
      info = data.frame(),
      genes = character()
    )
  )
  # Call select_genes with thresholds
  updated_object <- select_genes(
    angl = test_angl,
    zscore_mean_threshold = 10,
    zscore_sn_threshold = 10,
    max_n_genes = 3
  )

  # Verify that thresholds are adjusted and some genes are selected
  selected_genes <- updated_object@integration_genes$genes
  expect_true(length(selected_genes) > 0)
  expect_true(!any(is.na(selected_genes)))
})



