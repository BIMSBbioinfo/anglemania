## Don't forget to collect flippity metrics
## acess flippy status
flippalise <- function(l_processed) {
  ## make parallel TODO
  setNames(purrr::map(
    c("sharp", "blunt"),
    function(sset) {
      m <- l_processed[[paste0("x_", sset)]]
      df <- melt.to.df(m)
      df <- df[angle >= 3]
      df <- df[order(-angle)]
      df <- df[, edge := paste0(x, "_", y)]
      df$x <- NULL
      df$y <- NULL
      ##
      purrr::map2(
        l_processed$data_info,
        names(l_processed$data_info),
        function(data_info, ds) {
          df_ang <- data.table::fread(
            file = file.path(
              data_info$path_to_df_ang,
              data_info[[paste0("name_df_ang_", sset)]]
            )
          )
          df_ang <- df_ang[, edge := paste0(x, "_", y)]
          df_ang <- df_ang[, .(edge, flip.ref)]
          data.table::setnames(df_ang, old = "flip.ref", new = ds)
          df_ang_flp <- df_ang[df, on = "edge"]
          return(df_ang_flp)
        }
      ) %>% purrr::reduce(
        .,
        ~ data.table::merge.data.table(.x, .y, by = c("angle", "edge"))
      )
    }
  ) -> l_flp, c("sharp", "blunt"))
}
##
l_flp <- flippalise(l_processed)
view(l_flp$blunt[angle >= 8])
