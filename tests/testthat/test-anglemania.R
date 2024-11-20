load(system.file("extdata",
  "seurat_splatter_sim.RData",
  package = "anglemania"
))
se_raw <- se

test_that("anglemania works correctly with method pearson", {
  se <- se_raw
  anglemania_object <- create_anglemaniaObject(se, batch_key = "Batch")
  anglemania_object <- anglemania(anglemania_object, method = "pearson")

  # check that list_stats is not empty
  expect_true(length(anglemania_object@list_stats) > 0)

  # check that list_stats is a list
  expect_true(is.list(anglemania_object@list_stats))

  expect_snapshot(list_stats(anglemania_object)$mean_zscore[1:10, 1:10])

  # check if the first few elements from the list_stats
  # (zscore mean, zscore SD, signal-to-noise ratio) are the same
  expect_snapshot(str(list_stats(anglemania_object)))

  # check if the correct genes are extracted
  expect_snapshot(get_anglemania_genes(anglemania_object))
})

