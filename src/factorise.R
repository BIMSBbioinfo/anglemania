factorise <- function(x_mat, #nolint
                      extrema,
                      n_threads,
                      path_to_write_angles) {
  ## ---
  message(paste0("Input matrix with ", dim(x_mat)[1], " features"))
  if (is.null(rownames(x_mat))) {
    stop("Input matrix lacks rownames. Stopping")
  }
  ## --- extracting sample names and dimensionality for the record
  s_names <- colnames(x_mat)
  s_dims <- length(s_names)
  ## ---
  message("Computing cosine distances...")
  x_mat_ang <- extract.angles(x_mat, n_threads = n_threads)
  invisible(gc())
  ## ---
  message("Melting to data table...")
  x_df_ang <- melt.to.df(x_mat_ang)
  invisible(gc())
  ## ---
  message("Estimating critical angles...")
  l_angles <- estimate.critical.angles(x_df_ang, s_dims = s_dims, extrema)
  if (
    any(
      l_angles$critical_angles <= 0
    ) ||
      any(
        l_angles$critical_angles >= 180
      )
  ) {
    message(
      paste0("Angles at ",
        paste0(
               paste0(
                 l_angles$critical_angles, "%"
               ),
               collapse = " and "),
        ":\n",
        paste0(
          l_angles$critical_angles,
          collapse = " and "
        )
      )
    )
    stop("Extrema is too thin, consider taking a larger margin")
  }
  ## ---
  message("Assessing flippity...")
  l_flippity <- flippitise(x_mat, x_df_ang, l_angles)
  if (!is.na(path_to_write_angles) & is.character(path_to_write_angles)) {
    message("Writing tables with angles to compute integration metric...")
    l_df_ang_md5 <- write.angles(l_flippity, path_to_write_angles)
  }
  ## ---
  message("Factorising angle matrices...")
  x_mat_angfact_sharp <- ifelse(x_mat_ang <= l_angles$critical_angles[1], 1, 0)
  x_mat_angfact_sharp <- as(x_mat_angfact_sharp, "dgCMatrix")
  x_mat_angfact_blunt <- ifelse(x_mat_ang >= l_angles$critical_angles[2], 1, 0)
  x_mat_angfact_blunt <- as(x_mat_angfact_blunt, "dgCMatrix")
  invisible(gc())
  ##
  gib <- list(
    x_sharp = x_mat_angfact_sharp,
    x_blunt = x_mat_angfact_blunt,
    l_angles = l_angles,
    data_info = list(
      samples = s_names,
      path_to_df_ang = path_to_write_angles,
      name_df_ang_sharp = l_df_ang_md5[["sharp"]],
      name_df_ang_blunt = l_df_ang_md5[["blunt"]]
    )
  )
  return(gib)
}