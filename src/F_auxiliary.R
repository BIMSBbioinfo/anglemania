## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
extract.angles <- function(x_mat, #nolint
                           n_threads = 16) {
  x_mat <- as.matrix(x_mat)
  x_mat_ang <- round(
                     acos(
                       1 -
                         parallelDist::parDist(x_mat,
                                               method = "cosine",
                                               diag = FALSE,
                                               upper = FALSE,
                                               threads = 16)
                     ) / pi * 180,
                     2)
  x_mat_ang <- dist2mat(x_mat_ang, 256)
  return(x_mat_ang)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
melt.to.df <- function(x_mat_ang) { #nolint
  ##
  x_mat_ang[upper.tri(x_mat_ang, diag = TRUE)] <- NA
  x_df_ang <- data.table::as.data.table(x_mat_ang, keep.rownames = "x")
  x_df_ang <- data.table::melt(x_df_ang,
                               id.vars = "x",
                               variable.name = "y",
                               value.name    = "angle")
  x_df_ang <- na.omit(x_df_ang)
  x_df_ang <- x_df_ang[, y := as.character(y)]
  class(x_df_ang) <- "data.table"
  return(x_df_ang)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
aggregate.metrics <- function(l_processed) { #nolint
  purrr::map_dfr(l_processed$l_angles, ~ .x$statistics)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
get.plottab <- function(seu_obj, #nolint
                        var,
                        reduction) {
  cbind(
    seu_obj@meta.data %>% dplyr::select(orig.ident, dplyr::all_of(var)),
    seu_obj@reductions[[reduction]]@cell.embeddings
  )
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
min.max <- function(vec.num) { #nolint
  (vec.num - min(vec.num)) / (max(vec.num) - min(vec.num))
}