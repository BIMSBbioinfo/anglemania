test_that("Objects with fewer cells than min_mcells are removed", {
  result <- import_sl(sl_test_data, min_mcells = 6)
  expect_equal(length(result), 5)
})

test_that("SCTransform is applied if scaled data is not present", {
  result <- import_sl(sl_test_data, min_mcells = 6)
  all_true <- all(sapply(result, function(x) "scale.data" %in% SeuratObject::Layers(x)))
  expect_equal(all_true, TRUE)
})

# Tests for input validation
test_that("Error is thrown for list with fewer than two Seurat objects", {
  sl <- sl_test_data[1]
  expect_error(import_sl(sl, min_mcells = 6))
})

test_that("Error is thrown for non-Seurat objects in the list", {
  mock_seurat_list <- list("testtest", "asd")
  expect_error(import_sl(mock_seurat_list, min_mcells = 6))
})

test_that("Error is thrown for invalid min_mcells values", {
  expect_error(import_sl(sl_test_data, min_mcells = 5))
  expect_error(import_sl(sl_test_data, min_mcells = "invalid"))
})
