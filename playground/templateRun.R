## Template run for an integration task using ANGLEMANA
#### 24.08.2023
#### ab
######
source("./src/Header.R")
######
p_dat <- "/local/Projects/AAkalin_Neuroblastoma/Results/RDS/"
andat <- read_rds(
  file.path(p_dat,
    "AnnDataCsv_To_Seurat.be17bfa249d9d65a35d11e9fe4f1f6c7.rds"
  )
)
andat <- andat$dat
andat_l <- SplitObject(andat, split.by = "dataset")
workdat_l <- SplitObject(andat_l$Jansky, split.by = "sample_name")

# keep datasets with > 15 Scells
ds_keep <- purrr::map_dbl(names(workdat_l),
  ~ length(workdat_l[[.x]]$seacell_id)
) > 15
workdat_f_l <- workdat_l[ds_keep]
#  normalise & get gene expression
workdat_f_l <- lapply(workdat_f_l, FUN = NormalizeData)
workdat_f_l <- lapply(workdat_f_l, FUN = ScaleData)
l_mats <- purrr::map(workdat_f_l,
  function(pat) {
    pat@assays$RNA@scale.data
  }
)
# Process expression matrices
tic()
l_processed <- anglemanise(l_mats,
                           extrema = 0.001,
                           n_threads = 16,
                           path_to_write_angles = "./tmp")
toc()

# Extract features from angles ---- obsolete
tic()
l_features <- extract.integration.features(L_processed, cutoff = NULL)
toc()

# Integrate the dataset
tic()
seu_combined_test <- integrate.by.features(workdat_f_l, )
toc()