library(testthat)
library(data.table)
library(tidyr) # since approximate_angles uses tidyr
# Create a larger mock data table
set.seed(123) # for reproducibility
num_rows <- 100
x_df_ang <- data.table(
  x = paste0("gene", sample(1:50, num_rows, replace = TRUE)),
  y = paste0("gene", sample(1:50, num_rows, replace = TRUE)),
  angle = runif(num_rows, min = 0, max = 180) # random angles between 0 and 180 degrees
)
x_df_ang <- x_df_ang[x_df_ang$x != x_df_ang$x,]

test_that("approximate_angles handles inputs and outputs correctly", {
  # Create a mock data table
  
  # Test 1: Function returns error for non-positive quantile_split
  expect_error(approximate_angles(x_df_ang, quantile_split = 0))
  expect_error(approximate_angles(x_df_ang, quantile_split = -0.1))
  
  # Test 2: Function returns a list for valid inputs
  result <- approximate_angles(x_df_ang, quantile_split = 0.1)
  expect_type(result, "list")
  
  # Test 3: Check the structure of the returned list
  expect_true("critical_angles" %in% names(result))
  expect_true("angles_dist" %in% names(result))
  expect_true("angles_anylitical" %in% names(result))
  expect_true("statistics" %in% names(result))
  
  # Test 4: Validate the contents of 'angles_dist'
  expect_s3_class(result$angles_dist, "tbl_df")
  expect_true("angle" %in% colnames(result$angles_dist))
  expect_true("prob" %in% colnames(result$angles_dist))
  
  # Test 5: Ensure 'angles_dist' has expected values based on quantile_split
  expected_length <- 1 / 0.1 + 1  # number of intervals plus one
  expect_equal(nrow(result$angles_dist), expected_length)
  
  # Add more checks as necessary
})
