# load example seurat object
# pbmc_scmall does not have any "batches" from multiple experiments
#   but we'll just treat it like the "groups" column are the batches
se <- SeuratObject::pbmc_small
se_raw <- se

test_that("anglemania works correctly with method cosine", {
  se <- se_raw
  angl <- create_anglemania_object(se, batch_key = "groups")
  angl <- anglemania(angl, method = "cosine")

  # check that list_stats is not empty
  expect_true(length(angl@list_stats) > 0)

  # check that list_stats is a list
  expect_true(is.list(angl@list_stats))

  expect_snapshot(list_stats(angl)$mean_zscore[1:10, 1:10])

  # check if the first few elements from the list_stats
  # (zscore mean, zscore SD, signal-to-noise ratio) are the same
  expect_snapshot(list_stats(angl)$sn_zscore[1:10, 1:10])

  # check if the correct genes are extracted
  expect_snapshot(get_anglemania_genes(angl))
})

