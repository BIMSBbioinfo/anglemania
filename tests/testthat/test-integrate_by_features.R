load(system.file("extdata",
  "seurat_splatter_sim.RData",
  package = "anglemania"
))
se_raw <- se

test_that("integrate_by_features integrates Seurat objects correctly using selected features", {
  se <- se_raw
  anglemania_object <- create_anglemaniaObject(se, batch_key = "Batch")
  anglemania_object <- anglemania(anglemania_object, method = "pearson")

  # Integrate samples using selected features
  options(future.globals.maxSize = 5000 * 1024^2)
  suppressWarnings({
    se_integrated <- integrate_by_features(se, anglemania_object)
  })
  # Seurat gave too many unnecessary warnings.. Annoying in R-CMD-check
  
  # make snapshot of first few counts of integrated assay
  expect_snapshot(SeuratObject::LayerData(se_integrated)[1:10, 1:10])
})

