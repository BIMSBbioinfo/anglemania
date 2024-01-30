# test_that("input validation works correctly", {
#   expect_error(anglemanise(NULL), "seurat_list must be a list")
#   expect_error(anglemanise(list(), n_threads = "non-numeric"), "n_threads must be numeric")
#   expect_error(anglemanise(list(), n_threads = -1), "n_threads must be positive")
#   expect_error(anglemanise(list(), path_to_write_angles = 123), "path_to_write_angles must be a character string")
#   expect_error(anglemanise(list(), extrema = c(0.001, 0.002, 0.003)), "extrema must have length 1 or 2")
# })
# 
# 
# test_that("anglemanise returns expected structure", {
#   # Assuming you have a predefined 'seurat_list' for testing
#   result <- anglemanise(seurat_list, extrema = 0.001, n_threads = 2, path_to_write_angles = tempdir())
#   
#   expect_type(result, "list")
#   expect_true(all(c("x_sharp", "x_blunt", "l_angles", "data_info") %in% names(result)))
#   expect_true(is(result$x_sharp, "dgCMatrix")) #change to expect_type()
#   expect_true(is(result$x_blunt, "dgCMatrix"))
#   # Add more checks as necessary
# })
# 
# 
# test_that("anglemanise creates directory if non-existent", {
#   non_existent_path <- file.path(tempdir(), "new_dir")
#   expect_false(dir.exists(non_existent_path))
#   
#   anglemanise(seurat_list, extrema = 0.001, n_threads = 1, path_to_write_angles = non_existent_path)
#   expect_true(dir.exists(non_existent_path))
# })
# 
# 
# test_that("anglemanise writes files to specified directory", {
#   test_path <- file.path(tempdir(), "test_files")
#   dir.create(test_path)
#   
#   anglemanise(seurat_list, extrema = 0.001, n_threads = 1, path_to_write_angles = test_path)
#   
#   # Check if files are created. This will depend on how your function writes files
#   # Example:
#   # expect_true(length(list.files(test_path)) > 0)
# })
# 
