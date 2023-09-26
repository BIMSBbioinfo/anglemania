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

### Integrate the dataset by inferred features, iteratively reducing
### the featureset  by traversing the graph
l_seu_comb <- list()
l_seu_comb <- purrr::map(
    l_ifts,
    function(int_features) {
      integrate.by.features(workdat_f_l, int_features)
    }
)

### Now, pick the "best" feature set and integrate with inferred sample order
l_seu_comb$hclust_ord <- integrate.by.features(
  seurat_list = workdat_f_l[l_int_ordr$hclust$sample_order],
  features_to_integrate = l_ifts$all,
  int_order = l_int_ordr$hclust$int_matrix
)
l_seu_comb$graph_ord <- integrate.by.features(
  seurat_list = workdat_f_l[l_int_ordr$graph$sample_order],
  features_to_integrate = l_ifts$all,
  int_order = l_int_ordr$graph$int_matrix
)
l_seu_comb$graph_ord <- integrate.by.features(
  seurat_list = workdat_f_l[l_int_ordr$graph$sample_order],
  features_to_integrate = l_ifts$all,
  int_order = l_int_ordr$graph$int_matrix
)
l_seu_comb$graph_ord_rev <- integrate.by.features(
  seurat_list = workdat_f_l[rev(l_int_ordr$graph$sample_order)],
  features_to_integrate = l_ifts$all,
  int_order = l_int_ordr$graph$int_matrix
)
l_seu_comb$graph_ord_ved <- integrate.by.features(
  seurat_list = workdat_f_l[colnames(conm_summed)],
  features_to_integrate = l_ifts$all,
  int_order = l_int_ordr$graph$int_matrix
)



## First, show how reduced feature set affects integration
glist <- read_rds("./auxiliary/Glist_EMT.RDS")
purrr::map(
  l_seu_comb,
  ~ AddModuleScore(.x, glist, name = 'module', seed = 42)
) -> l_seu_comb

get.plottab <- function(seu_obj, var, reduction) {
  cbind(
    seu_obj@meta.data %>% dplyr::select(orig.ident, all_of(var)),
    seu_obj@reductions[[reduction]]@cell.embeddings
  )
}
purrr::map(
  l_seu_comb,
  ~ get.plottab(.x, "module1", "umap")
) -> l_seu_ptabs
l_seu_plots <- purrr::map2(
  l_seu_ptabs,
  names(l_seu_ptabs),
  ~ ggplot(data = .x, aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(color = module1), size = 2.5, alpha = 0.75) + 
  scale_colour_gradientn(name = "EMT score", colours = viridis::viridis(100)) + 
  labs(title = .y) +
  theme_minimal()
)
l_seu_plots$graph_ord_ved
