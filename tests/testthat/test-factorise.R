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
  permuted_values <- permuted_correlation[upper.tri(
    permuted_correlation,
    diag = FALSE
  )]
  mean_permuted <- mean(permuted_values)
  sd_permuted <- sd(permuted_values)

  # Step 7: Transform original correlation matrix into z-scores
  expected_z_scores <- (original_correlation - mean_permuted) / sd_permuted
  diag(expected_z_scores) <- NA
  # Compare the result from factorise with the expected z-score matrix
  expect_equal(result_matrix, expected_z_scores)
})


test_that("factorise computes z-score-transformed 
angle matrix correctly with method 'diem'", {
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
  
  # Run factorise with method "diem" and a fixed seed
  result_fbm <- factorise(x_mat, method = "diem", seed = 1)
  
  # Extract the result as a regular matrix
  result_matrix <- result_fbm[]
  
  # Manually compute the expected result
  # Step 1: Permute the matrix column-wise with the same seed
  set.seed(1)
  permuted_data <- apply(mat, 2, sample)
  
  # Step 2: Normalize and scale both original and
  # permuted data using Seurat::NormalizeData


  normalized <- Seurat::NormalizeData(
    mat,
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )
  permuted_normalized <- Seurat::NormalizeData(
    permuted_data,
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )

  
  # Step 3: Compute cosine correlation matrices
  # Transpose the normalized data to have genes as columns for cor()
  original_transposed <- t(normalized)
  permuted_transposed <- t(permuted_normalized)
  
  original_correlation <- cor(
    original_transposed,
    use = "pairwise.complete.obs"
  )
  permuted_correlation <- cor(
    permuted_transposed,
    use = "pairwise.complete.obs"
  )
  
  # Step 4: Set diagonal to NA
  diag(original_correlation) <- NA
  diag(permuted_correlation) <- NA
  
  # Step 5: Compute Euclidean distances for DIEM: sqrt(2 * (1 - r))
  diem_correlation <- sqrt(2 * (1 - original_correlation))
  
  # Step 6: Compute dstat from the permuted correlation matrix
  # Assuming get_dstat is accessible and computes mean, var, sd
  # Step 7: Compute scale_factor
  scale_factor <- (min(diem_correlation, na.rm = TRUE) -
    max(diem_correlation, na.rm = TRUE)) /
    var(diem_correlation[upper.tri(diem_correlation, diag = FALSE)])
  
  # Step 8: Scale the DIEM distances
  scaled_diem <- scale_factor * (diem_correlation -
    mean(diem_correlation[upper.tri(diem_correlation, diag = FALSE)]))
  
  # Step 9: Transform to z-scores
  expected_z_scores <- (scaled_diem -
    mean(diem_correlation[upper.tri(diem_correlation, diag = FALSE)])) /
    sd(diem_correlation[upper.tri(diem_correlation, diag = FALSE)])
  
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
  permuted_values <- permuted_correlation[upper.tri(permuted_correlation, diag = FALSE)]
  mean_permuted <- mean(permuted_values)
  sd_permuted <- sd(permuted_values)

  # Step 7: Transform original correlation matrix into z-scores
  expected_z_scores <- (original_correlation - mean_permuted) / sd_permuted
  diag(expected_z_scores) <- NA
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
