## Report on the assessment of integration orders on the example
## of Jansky and Dong datasets
#### 26.09.2023
#### ab
######

###### INFO
## Datasets will be integrated using direct/reversed orders derived from
## connectivity graph and dendrogram based on the same similarity metric.

## To have some notion of whats going on, EMT and Cell-Cycle will
## be scored independently across meta-cells.


source("./src/Header.R")
######
for (ds in c("Bedoya", "Dong", "Jansky-Adrenal", "Jansky")) {
  sl <- import.sl(path_dat = paste0("data/SL_", ds, ".rds"))
  l_processed <- anglemanise(
    sl,
    extrema = 0.005,
    n_threads = 16,
    path_to_write_angles = "./tmp"
  )
  readr::write_rds(
    l_processed,
    file = paste0("results/", ds, "_processed_0.01.rds")
  )
}

## prepare lists of individual features per intersected ds
featurised <- featurise(l_processed)
plot.unqiue.fts(featurised$unique_features, n_added)




