test_that("extract_angles extracts cosine similarities correctly", {
  mat <- matrix(
    c(0, 0, 0, 4, 2,
      1, 2, 3, 4, 5,
      5, 4, 3, 2, 1,
      0, 0, 0, 0, 0),
    nrow = 4,
    ncol = 5,
    byrow = TRUE
  )
  
  # Convert to FBM
  x_mat <- bigstatsr::FBM(nrow = nrow(mat), ncol = ncol(mat), init = mat)
  
  # Run extract_angles with method "cosine"
  result_fbm <- extract_angles(x_mat, method = "cosine")
  
  # Extract the result as a regular matrix
  result_matrix <- result_fbm[]
  
  # Manually compute the expected result
  log_normalized_data <- Seurat::NormalizeData(mat)
  # normalized_data <- t(t(mat) / colSums(mat) * 10000)
  # log_normalized_data <- log1p(normalized_data)
  
  transposed_data <- t(log_normalized_data)
  suppressWarnings({
    expected_correlation <- cor(transposed_data, use = "pairwise.complete.obs")
  })
  
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
      0, 0, 0, 0, 0
    ),
    nrow = 4,
    ncol = 5,
    byrow = TRUE
  )

  # Convert to FBM
  x_mat <- bigstatsr::FBM(nrow = nrow(mat), ncol = ncol(mat), init = mat)

  # Run extract_angles with method "cosine"
  result_fbm <- extract_angles(x_mat, method = "spearman")

  # Extract the result as a regular matrix
  result_matrix <- result_fbm[]

  # Manually compute the expected result
  log_normalized_data <- Seurat::NormalizeData(mat)
  # normalized_data <- t(t(mat) / colSums(mat) * 10000)
  # log_normalized_data <- log1p(normalized_data)

  transposed_data <- t(log_normalized_data)
  suppressWarnings({
    expected_correlation <- cor(
      transposed_data,
      method = "spearman"
    )
  }) # warning thrown because of standard deviation equal to zero but
  # we allow for this because t

  diag(expected_correlation) <- NA

  expect_s4_class(result_fbm, "FBM")
  # Compare the result
  expect_equal(result_matrix, expected_correlation)
})