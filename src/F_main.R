## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
write.angles <- function(l_flippity, #nolint
                         path_to_write_angles) {
  ##
  l_df_ang_md5 <- list()
  for (i in names(l_flippity)) {
    md5 <- digest::digest(l_flippity[[i]], "md5")
    data.table::fwrite(
      x = l_flippity[[i]],
      file = file.path(path_to_write_angles, md5), 
      sep = "\t"
    )
    l_df_ang_md5[[i]] <- md5
  }
  return(l_df_ang_md5)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## Get features from cumulative factang matrices assess distributions
## Currently Under BIIIG question
extract.integration.features <- function(l_cumfact, #nolint
                                         cutoff = NULL) {
  purrr::map(l_cumfact[1:2],
             function(cumfact_mat) {
               cumfact_dt <- as.data.table(cumfact_mat, keep.rownames = "x")
               cumfact_dt_melt <- melt(cumfact_dt, 
                                       id.vars = "x",
                                       variable.name = "y",
                                       value.name    = "angle.fact"
               )
               ## Clean the table
               cumfact_dt_melt <- cumfact_dt_melt[!is.na(angle.fact) & (x != y)]
               rm(cumfact_dt); gc()
               ## this distribution likely follows poisson ???; get some basis to select a cutoff
               if (is.null(cutoff)) {
                 cumfact_dt_melt %>% 
                   tabyl(angle.fact) %>% 
                   as_tibble() -> fact_stats
                 glm(n ~ angle.fact, data = fact_stats, family = quasi(variance = "mu", link = "log")) -> glmodel.poisson
                 fact_stats %>%
                   mutate(n_pred_poisson = predict(glmodel.poisson, 
                                                   newdata = tibble(angle.fact = fact_stats$angle.fact)) %>% 
                            exp() %>%
                            unname() %>% 
                            as.integer()
                   ) -> fact_stats_pred
                 fact_stats_pred %>% 
                   transmute(diff = log(n/n_pred_poisson),
                             diff = ifelse(diff == Inf, 0, diff)) %>% 
                   pull(diff) -> v_diff
                 which(v_diff == max(v_diff)) - 1 -> cutoff
               } # first -1 to get to 0-based counting; 
               ##
               cumfact_dt_melt_filt <- cumfact_dt_melt[ angle.fact >= cutoff ]
               rm(cumfact_dt_melt); gc()
               features <- intersect(cumfact_dt_melt_filt$x, cumfact_dt_melt_filt$y)
               ##
               return(features)
             }
  ) -> L_features
  ##
  names(L_features) <- names(L_cumfact[1:2])
  #L_features$x_blunt <- L_features$x_blunt[ !L_features$x_blunt %in% purrr::reduce(L_features, intersect) ]
  return(L_features)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## for now, just use CCA and Vedran's settings with Seurat
integrate.by.features <- function(seurat_list, features_to_integrate, int_order = NULL) {
  anchors = FindIntegrationAnchors(
    object.list = seurat_list, 
    anchor.features = features_to_integrate, 
    verbose = FALSE,
    dims = 1:10,
    k.filter  = 10, 
    k.anchor  = 10, 
    k.score   = 10,
    reduction = "cca"
  )
  features_intersect = Reduce(function(x,y) intersect(x,y), lapply(seurat_list, rownames))
  seurat_combined = IntegrateData(
    anchorset = anchors,
    k.weight = 10, 
    features.to.integrate = features_intersect,
    sample.tree = int_order,
    verbose = T
  )
  DefaultAssay(seurat_combined) = "integrated"
  seurat_combined <- ScaleData(seurat_combined)
  seurat_combined <- RunPCA(seurat_combined, npcs = 30)
  seurat_combined <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:30)
  return(seurat_combined)
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
extract.integration.statistics <- function(Seurat_integrated, L_factmats) {
  ## Get angles for the integrated dataset
  X_int <- Seurat_integrated@assays$integrated@scale.data
  #X_intrg_ang <- extract.angles(X_intrg)
  #X_intrg_ang_df <- melt.to.df(X_intrg_ang)
  #X_intrg_ang_df <- X_intrg_ang_df[ x != y ]
  ##
  purrr::map2_dfr(L_processed$data_info,
                  names(L_processed$data_info),
                  function(dataset, data_name) {
                    X_int_ssmp <- X_int[, dataset$samples]
                    X_int_ssmp_ang <- extract.angles(X_int_ssmp)
                    X_int_ssmp_ang_df <- melt.to.df(X_int_ssmp_ang)
                    X_int_ssmp_ang_df <- X_int_ssmp_ang_df[ x != y ]
                    rm(X_int_ssmp, X_int_ssmp_ang); gc()
                    ##
                    sharp.mrgd <- collect.written.angles(X_int_ssmp_ang_df, dataset, data_name, "sharp")
                    blunt.mrgd <- collect.written.angles(X_int_ssmp_ang_df, dataset, data_name, "blunt")
                    ##
                    out <- rbind(sharp.mrgd, blunt.mrgd)
                    return(out)
                  }) -> angdif_df
  ##
  return(angdif_df)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
collect.written.angles <- function(X_postint_ang_df, dataset, data_name, ang_type) {
  ## 
  tmp <- data.table::fread(file = file.path(dataset$path_to_df_ang, 
                                            dataset[[paste0("name_df_ang_", ang_type)]]
                                            )
                           )
  tmp[X_postint_ang_df, on = c("x", "y"), nomatch = 0] -> tmp.mrgd
  colnames(tmp.mrgd) <- c("x", "y", "angle_prei", "angle_posti")
  ## summarised
  tmp.mrgd[, .(angle.mean = mean(angle_prei), angle.sd  = sd(angle_prei))] -> prei
  tmp.mrgd[, .(angle.mean = mean(angle_posti), angle.sd  = sd(angle_posti))] -> posti
  rbind(
    prei[, c("angle.source", "angle_prei.type", "data.source") := .("prei", ang_type, data_name)],
    posti[, c("angle.source", "angle_prei.type", "data.source") := .("posti", ang_type, data_name)]
  ) -> tmp.mrgd
  return(tmp.mrgd)
}