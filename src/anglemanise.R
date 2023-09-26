anglemanise <- function(seurat_list, #nolint
                        extrema = 0.001,
                        n_threads = 16,
                        path_to_write_angles = ".") {
  ##
  if (is.numeric(extrema)) {
    if (length(extrema) > 1) {
      if (length(extrema) > 2) {
        stop(paste0("extrema should be numeric of length either 1 or 2"))
      }
      extrema <- extrema
    } else {
      extrema <- c(extrema, 1 - extrema)
    }
  } else {
    stop("extrema should be numeric of length either 1 or 2")
  }
  #-------------
  message(
    paste0("Processing ",
           length(list_x_mats),
           " expression matrices...")
    )
  ##
  if (!dir.exists(path_to_write_angles)) {
    dir.create(path_to_write_angles, recursive = TRUE)
  }
  ## - extract expression matrices
  list_x_mats <- purrr::map(
    seurat_list,
    function(ss) {
      ss@assays$RNA@scale.data
    }
  )
  ##
  if (is.null(names(list_x_mats))) {
    names(list_x_mats) <- paste0("X", seq_along(list_x_mats))
  }
  ##
  g_dims <- dim(list_x_mats[[1]])[1]
  l_added <- list(
    x_sharp = Matrix::Matrix(0, nrow = g_dims, ncol = g_dims, sparse = TRUE),
    x_blunt = Matrix::Matrix(0, nrow = g_dims, ncol = g_dims, sparse = TRUE),
    data_info = list()
  )
  ##
  p <- progressr::progressor(along = list_x_mats)
  for (mat_ind in seq_along(list_x_mats)) {
    ##
    x_name <- names(list_x_mats[mat_ind])
    ##
    message(paste0("Starting matrix ", x_name))
    Sys.sleep(1)
    ##
    p(message = sprintf("Processing %s", names(list_x_mats)[mat_ind]))
    x <- factorise(list_x_mats[[x_name]],
                   extrema,
                   n_threads,
                   path_to_write_angles)
    ##
    message(paste0("Adding ", x_name))
    l_added[["x_sharp"]] <- l_added[["x_sharp"]] + x[["x_sharp"]]
    l_added[["x_blunt"]] <- l_added[["x_blunt"]] + x[["x_blunt"]]
    l_added[["l_angles"]][[x_name]] <- x[["l_angles"]]
    l_added[["data_info"]][[x_name]] <- x[["data_info"]]
    ##
    invisible(gc())
  }
  ## Implement parallel --- heavy memory load, but potential big speed-up
  ### using future



  ##
  return(l_added)
}