test_that("factorise computes z-score-transformed angle 
matrix correctly with method 'cosine'", {
  # Create a small test matrix with known values (genes x cells)
  # Let's assume 4 genes and 5 cells
  mat <- matrix(
    c(
      5, 3, 0, 0, 2,
      0, 0, 0, 3, 4,
      2, 1, 3, 4, 0,
      0, 0, 1, 0, 0
    ),
    nrow = 4,
    ncol = 5,
    byrow = TRUE
  )

  # Convert to FBM
  set.seed(1)
  x_mat <- bigstatsr::FBM(nrow = nrow(mat), ncol = ncol(mat), init = mat)

  # Run factorise with method "cosine" and a fixed seed
  result_fbm <- factorise(x_mat, method = "cosine", seed = 1)

  # Extract the result as a regular matrix
  result_matrix <- result_fbm[]

  # Manually compute the expected result

  # Step 1: Permute the matrix column-wise
  set.seed(1)
  permuted_data <- apply(mat, 2, sample)

  # Step 2: Normalize and scale both original and permuted data
  log_normalized <- Seurat::NormalizeData(mat)
  permuted_log_normalized <- Seurat::NormalizeData(permuted_data)


  # Step 3: Compute cosine correlation matrices
  original_correlation <- cor(
    t(log_normalized),
    use = "pairwise.complete.obs"
  )
  permuted_correlation <- cor(
    t(permuted_log_normalized),
    use = "pairwise.complete.obs"
  )


  # Step 4: Compute mean and standard deviation
  # from the permuted correlation matrix
  # Exclude NA values
  diag(permuted_correlation) <- NA
  mean_permuted <- Matrix::colMeans(permuted_correlation, na.rm = TRUE)
  sd_permuted <- matrixStats::colSds(permuted_correlation, na.rm = TRUE)

  # Step 7: Transform original correlation matrix into z-scores
  expected_z_scores <- (original_correlation - mean_permuted) / sd_permuted
  diag(expected_z_scores) <- 0
  # Compare the result from factorise with the expected z-score matrix
  expect_equal(result_matrix, expected_z_scores)
})


test_that("factorise computes z-score-transformed angle
matrix correctly with method 'spearman'", {
  # Create a small test matrix with known values (genes x cells)
  # Let's assume 4 genes and 5 cells
  mat <- matrix(
    c(
      5, 3, 0, 0, 2,
      0, 0, 0, 3, 4,
      2, 1, 3, 4, 0,
      0, 0, 1, 0, 0
    ),
    nrow = 4,
    ncol = 5,
    byrow = TRUE
  )

  # Convert to FBM
  x_mat <- bigstatsr::FBM(nrow = nrow(mat), ncol = ncol(mat), init = mat)

  # Run factorise with method "cosine" and a fixed seed
  result_fbm <- factorise(x_mat, method = "spearman", seed = 1)

  # Extract the result as a regular matrix
  result_matrix <- result_fbm[]

  # Manually compute the expected result

  # Step 1: Permute the matrix column-wise
  set.seed(1)
  permuted_data <- apply(mat, 2, sample)

  # Step 2: Normalize and scale both original and permuted data
  log_normalized <- Seurat::NormalizeData(mat)
  permuted_log_normalized <- Seurat::NormalizeData(permuted_data)


  # Step 3: Compute cosine correlation matrices
  original_correlation <- cor(
    t(log_normalized),
    use = "pairwise.complete.obs",
    method = "spearman"
  )
  permuted_correlation <- cor(
    t(permuted_log_normalized),
    use = "pairwise.complete.obs",
    method = "spearman"
  )


  # Step 4: Compute mean and standard deviation
  # from the permuted correlation matrix
  # Exclude NA values
  diag(permuted_correlation) <- NA
  mean_permuted <- Matrix::colMeans(permuted_correlation, na.rm = TRUE)
  sd_permuted <- matrixStats::colSds(permuted_correlation, na.rm = TRUE)

  # Step 7: Transform original correlation matrix into z-scores
  expected_z_scores <- (original_correlation - mean_permuted) / sd_permuted
  diag(expected_z_scores) <- 0
  # Compare the result from factorise with the expected z-score matrix
  expect_equal(result_matrix, expected_z_scores)
})


# Edge cases
test_that("factorise handles empty matrices correctly", {
  # Empty matrix with zero rows
  empty_mat_rows <- matrix(nrow = 0, ncol = 5)
  x_mat_rows <- bigstatsr::FBM(nrow = 0, ncol = 5)
  
  expect_error(factorise(x_mat_rows, method = "cosine"))
  
  # Empty matrix with zero columns
  empty_mat_cols <- matrix(nrow = 5, ncol = 0)
  x_mat_cols <- bigstatsr::FBM(nrow = 5, ncol = 0)
  
  expect_error(factorise(x_mat_cols, method = "cosine"))
})

test_that("factorise really computes the gene-gene correlation matrix -
check dimensions", {
  # Create a small test matrix with known values (genes x cells)
  # Let's assume 4 genes and 5 cells
  mat <- matrix(
    c(
      5, 3, 0, 0,
      0, 0, 0, 3,
      2, 1, 3, 4,
      0, 0, 1, 0,
      1, 2, 1, 2,
      3, 4, 3, 4
    ),
    nrow = 6, # 6 genes
    ncol = 4, # 4 cells
    byrow = TRUE
  )

  mat <- bigstatsr::FBM(nrow = nrow(mat), ncol = ncol(mat), init = mat)


  # Run factorise with method "cosine" and a fixed seed
  result_fbm <- factorise(mat, method = "cosine", seed = 1)
  
  # check dims
  expect_equal(nrow(result_fbm), nrow(mat))
  expect_equal(ncol(result_fbm), nrow(mat))
})
