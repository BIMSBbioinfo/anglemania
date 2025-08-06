sce_raw <- sce_example()
test_that("anglemania works correctly with method cosine", {
  library(S4Vectors)
  sce <- sce_raw
  sce <- anglemania(sce, batch_key = "batch")

  # check that list_stats is not empty
  expect_true(length(metadata(sce)$anglemania$list_stats) > 0)

  # check that list_stats is a list
  expect_true(is.list(metadata(sce)$anglemania$list_stats))

  expect_snapshot(metadata(sce)$anglemania$list_stats$mean_zscore[1:10, 1:10])

  # check if the first few elements from the list_stats
  # (zscore mean, zscore SD, signal-to-noise ratio) are the same
  expect_snapshot(metadata(sce)$anglemania$list_stats$sn_zscore[1:10, 1:10])

  # check if the correct genes are extracted
  expect_snapshot(get_anglemania_genes(sce))
})
