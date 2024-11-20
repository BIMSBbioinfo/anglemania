test_that("extract_angles extracts pearson correlations correctly", {
  mat <- matrix(
    c(0, 0, 0, 4, 2,
      1, 2, 3, 4, 5,
      5, 4, 3, 2, 1,
      0, 1, 0, 1, 0),
    nrow = 4,
    ncol = 5,
    byrow = TRUE
  )
  
  # Convert to FBM
  x_mat <- bigstatsr::FBM(nrow = nrow(mat), ncol = ncol(mat), init = mat)
  
  # Run extract_angles with method "pearson"
  result_fbm <- extract_angles(x_mat, method = "pearson")
  
  # Extract the result as a regular matrix
  result_matrix <- result_fbm[]
  
  # Manually compute the expected result
  log_normalized_data <- Seurat::NormalizeData(mat)
  # normalized_data <- t(t(mat) / colSums(mat) * 10000)
  # log_normalized_data <- log1p(normalized_data)
  
  transposed_data <- t(log_normalized_data)
  expected_correlation <- cor(transposed_data, use = "pairwise.complete.obs")
  
  diag(expected_correlation) <- NA
  
  expect_s4_class(result_fbm, "FBM")
  # Compare the result
  expect_equal(result_matrix, expected_correlation)

})

test_that("extract_angles extracts spearman correlations correctly", {
  mat <- matrix(
    c(
      0, 0, 0, 4, 2,
      1, 2, 3, 4, 5,
      5, 4, 3, 2, 1,
      0, 1, 0, 1, 0
    ),
    nrow = 4,
    ncol = 5,
    byrow = TRUE
  )

  # Convert to FBM
  x_mat <- bigstatsr::FBM(nrow = nrow(mat), ncol = ncol(mat), init = mat)

  # Run extract_angles with method "pearson"
  result_fbm <- extract_angles(x_mat, method = "spearman")

  # Extract the result as a regular matrix
  result_matrix <- result_fbm[]

  # Manually compute the expected result
  log_normalized_data <- Seurat::NormalizeData(mat)
  # normalized_data <- t(t(mat) / colSums(mat) * 10000)
  # log_normalized_data <- log1p(normalized_data)

  transposed_data <- t(log_normalized_data)
  expected_correlation <- cor(
    transposed_data, 
    use = "pairwise.complete.obs",
    method = "spearman"
  )

  diag(expected_correlation) <- NA

  expect_s4_class(result_fbm, "FBM")
  # Compare the result
  expect_equal(result_matrix, expected_correlation)
})
