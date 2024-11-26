# load example seurat object
# pbmc_scmall does not have any "batches" from multiple experiments
#   but we'll just treat it like the "groups" column are the batches
se <- SeuratObject::pbmc_small
se_raw <- se

test_that("integrate_by_features integrates Seurat objects correctly using selected features", {
  se <- se_raw
  angl <- create_anglemania_object(se, batch_key = "groups")
  angl <- anglemania(angl, method = "pearson")

  # Integrate samples using selected features
  options(future.globals.maxSize = 5000 * 1024^2)
  suppressWarnings({
    se_integrated <- integrate_by_features(se, angl)
  })
  # Seurat gave too many unnecessary warnings.. Annoying in R-CMD-check
  
  # make snapshot of first few counts of integrated assay
  expect_snapshot(SeuratObject::LayerData(se_integrated)[1:10, 1:10])
})

