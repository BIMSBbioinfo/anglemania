context("anglemanise function tests")

test_that("Input validation works correctly", {
  # Invalid n_threads
  expect_error(anglemanise(sl_scaled_test_data, n_threads = -1))
  
  # Invalid path_to_write_angles
  expect_error(anglemanise(sl_scaled_test_data, path_to_write_angles = 123))
})


result <- anglemanise(sl_scaled_test_data)

test_that("Function handles valid inputs correctly", {
  # Verify the structure of the result
  expect_type(result, list)
  expect_equal(length(result), 4)
  expect_true("x_sharp" %in% names(result))
  expect_true("x_blunt" %in% names(result))
  expect_true("data_info" %in% names(result))
  expect_true("l_angles" %in% names(result))
  
})

test_that("Sharp and Blunt matrices are consistent", {
  sharp_data <- matrix(c(
    5, 0, 0, 1, 0, 1, 0, 1, 0, 1,
    0, 5, 1, 0, 0, 0, 0, 0, 0, 2,
    0, 1, 5, 2, 0, 2, 0, 0, 0, 3,
    1, 0, 2, 5, 1, 1, 1, 0, 1, 1,
    0, 0, 0, 1, 5, 0, 1, 1, 0, 1,
    1, 0, 2, 1, 0, 5, 1, 0, 0, 1,
    0, 0, 0, 1, 1, 1, 5, 1, 0, 0,
    1, 0, 0, 0, 1, 0, 1, 5, 0, 1,
    0, 0, 0, 1, 0, 0, 0, 0, 5, 0,
    1, 2, 3, 1, 1, 1, 0, 1, 0, 5),
    nrow = 10, byrow = TRUE)
  sharp_sparse_matrix <- as(data, "dgCMatrix")
  sharp_sparse_matrix@Dimnames[[2]] <- c("AL669831.5", "NOC2L", "HES4", "ISG15",
                                         "SDF4", "UBE2J2", "ACAP3", "INTS11",
                                         "AURKAIP1", "CCNL2")
  
  blunt_data <- matrix(c(
    0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 2, 0, 1, 1, 1, 0,
    1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 2, 0, 0, 0, 0, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
    nrow = 10, byrow = TRUE)
  blunt_sparse_matrix <- as(data, "dgCMatrix")
  blunt_sparse_matrix@Dimnames[[2]] <- c("AL669831.5", "NOC2L", "HES4", "ISG15",
                                         "SDF4", "UBE2J2", "ACAP3", "INTS11", 
                                         "AURKAIP1", "CCNL2")
  
  first_100_elements_sharp <- result$x_sharp[1:10,1:10]
  first_100_elements_blunt <- result$x_blunt[1:10,1:10]
  
  expect_snapshot(
    waldo::compare(first_100_elements_sharp, sharp_sparse_matrix)
  )
  expect_snapshot(
    waldo::compare(first_100_elements_blunt, sharp_sparse_matrix_blunt)
  )
  })


test_that("Check anglemanised $angles_dist", {
  # Verify the structure of the result
  expect_type(result, list)
  expect_equal(length(result), 4)
  expect_true("x_sharp" %in% names(result))
  expect_true("x_blunt" %in% names(result))
  expect_true("data_info" %in% names(result))
  expect_true("l_angles" %in% names(result))
  
})


test_that("Function handles edge cases correctly", {
  # Edge cases like empty Seurat objects, objects with no scaled data, etc.
  mock_seurat_list <-  list("testtest", "asd")
  expect_error(anglemanise(mock_seurat_list))
})

