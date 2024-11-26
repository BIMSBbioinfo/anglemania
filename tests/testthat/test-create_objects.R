# load example seurat object
# pbmc_scmall does not have any "batches" from multiple experiments
#   but we'll just treat it like the "groups" column are the batches
se <- SeuratObject::pbmc_small
se_raw <- se

# Should check if all the slots are correctly
## constructed and in the end check the entire angl object

## unique batch key "batch" is added
test_that("unique batch key is added", {
  se <- se_raw
  batch_key <- "groups"
  dataset_key <- "Dataset"
  dataset_col <- data.frame(
    Dataset = rep(c("Dataset1", "Dataset2"), ncol(se) / 2),
    row.names = colnames(se)
  )
  se <- Seurat::AddMetaData(se, dataset_col)

  se <- add_unique_batch_key(se,
    batch_key = batch_key,
    dataset_key = dataset_key
  )
  meta <- se[[]]

  expect_true("batch" %in% colnames(meta))

})


## matrix_list
test_that("matrix list is correctly constructed", {
  library(Seurat)
  batch_key <- "groups"
  matrices <- Seurat::SplitObject(se, split.by = batch_key) |>
    lapply(SeuratObject::LayerData, layer = "counts")
  matrices <- lapply(matrices, function(x) {
    x <- sparse_to_fbm(x)
    x
  })

  angl <- create_anglemania_object(se, batch_key = batch_key)
  expect_true(length(matrix_list(angl)) == 2)
  expect_s4_class(matrix_list(angl)[[1]], "FBM")
})

test_that("data_info and weights correctly constructed 
in case of no dataset_key", {
  se <- se_raw
  batch_key <- "groups"
  data_info_df <- se[[]] %>%
    dplyr::mutate(batch = groups) %>%
    dplyr::select(batch, all_of(batch_key)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(weight = 1 / nrow(.))
  
  # check data_info
  angl <- create_anglemania_object(se, batch_key = batch_key)
  data_info <- data_info(angl)
  expect_equal(data_info, data_info_df)

  # check weights
  weights <- angl_weights(angl)
  real_weights <- data_info_df$weight
  names(real_weights) <- data_info_df$batch
  expect_equal(weights, real_weights)

  # check dataset_key
  expect_s4_class(dataset_key(angl), NA)
})

test_that(
  "data_info and weights correctly constructed in case of 
  existing dataset_key",
  {
    se <- se_raw
    batch_key <- "groups"
    dataset_key <- "Dataset"
    dataset_col <- data.frame(
      Dataset = rep(c("Dataset1", "Dataset2"), ncol(se) / 2),
      row.names = colnames(se)
    )
    se <- SeuratObject::AddMetaData(se, dataset_col)

    se <- add_unique_batch_key(se,
      batch_key = batch_key,
      dataset_key = dataset_key
    )
    data_info_df <- se[[]] %>%
      dplyr::select(
        batch,
        all_of(c(dataset_key, batch_key))
      ) %>%
      dplyr::distinct() %>%
      dplyr::group_by(across(all_of(dataset_key))) %>%
      dplyr::add_count(across(all_of(dataset_key)), name = "n_samples") %>%
      dplyr::mutate(
        weight = 1 / n_samples / dplyr::n_groups(.)
      )
    # check data_info
    angl <- create_anglemania_object(
      se,
      batch_key = batch_key,
      dataset_key = dataset_key
    )
    data_info <- data_info(angl)
    expect_equal(data_info, data_info_df)
  }
)

test_that("anglemania_object is correct in case of incorrect dataset_key", {
  library(Seurat)
  library(dplyr)
  se <- se_raw
  batch_key <- "groups"
  dataset_key <- "false"
  expect_error(create_anglemania_object(se,
    batch_key = batch_key,
    dataset_key = dataset_key
  ))
  dataset_key <- 1
  expect_error(create_anglemania_object(se,
    batch_key = batch_key,
    dataset_key = dataset_key
  ))
})
