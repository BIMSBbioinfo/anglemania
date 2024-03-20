#' Factorise angle matrices
#'
#' @description
#' Generates two a one-hot encoded square matrixces
#' recording the presence/absence (1/0) of sharp and blunt
#' critical angles.
#'
#' @details
#' *factorise* extracts angles between genes, estimates
#' critical angles by approximation of the angle distribution,
#' and records angles passing the critical threshold into
#' a sparse matrix.
#'
#' @param x_mat Matrix. Contains normalised and scaled gene
#'   expression.
#' @param extrema double. Fraction of the angles
#'   to be cut from both sides of an approximated angle
#'   distribution.
#' @param n_threads integer. Number of cores to use for
#'   the computation of the cosine distances.
#' @param path_to_write_angles string. A path, under which
#'   a temporary directory will be created to record values
#'   of critical angles for each sample.
#' @return list. First two elements are sparse matrices
#'   containing the significant sharp and blunt angles
#'   between genes. Third and fourth elements are lists
#'   with angles statistics and paths to the values of
#'   critical mangles.
#' @export factorise
factorise <- function(x_mat, #nolint
                      extrema,
                      n_threads,
                      path_to_write_angles) {
  if (is.null(rownames(x_mat))) {
    stop("Input matrix lacks rownames. Stopping")
  }
  s_names <- colnames(x_mat)
  s_dims <- length(s_names)
  message("Computing cosine distances...")
  x_mat_ang <- extract_angles(x_mat, n_threads = n_threads)
  invisible(gc())
  # NOTE: For now, don't use melt_to_df, filter_angles, and write_anlges because it's slow
  # NOTE: they belong to the assemble cons nodes etc. code and bascially store the angle data.tables for each sample, 
  # then later, the data.tables  are fioltered by the conserved gene pairs of the processed x_mat_ang. Then it is basically 
  # computes which data.tables/matrices are the most similar, i.e., how many blunt/sharp gene pairs each sample pair share
  # message("Melting to data table...")
  # x_df_ang <- melt_to_df(x_mat_ang)
  invisible(gc())
  message("Estimating critical angles...")
  l_angles <- estimate_critical_angles(x_mat_ang, s_dims = s_dims, extrema)
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
  message("Filtering angles...")
  ##### DEPRECATE WRITING ANGLES FOR NOW AND SEE HOW IT GOES?
  # FIXME: taking too long.
  # COMMENT: an idea was to instead of storing and tracing back which samples are the closest, one could also use a binary matrix 
  # COMMENT: https://chat.openai.com/share/34bc83e7-ad3a-49aa-a7cc-a08fb4e89861 
  # COMMENT: (based on bytes instead of 0 and 1) to record which samples contributed to which element in the x_mat_ang.
  # l_flippity <- filter_angles(x_df_ang, l_angles)
  # if (!is.na(path_to_write_angles) & is.character(path_to_write_angles)) {
  #   message("Writing tables with angles to compute integration metric...")
  #   l_df_ang_md5 <- write_angles(l_flippity, path_to_write_angles)
  # }
  message("Factorising angle matrices...")
  x_mat_angfact_sharp <- ((x_mat_ang <= l_angles$critical_angles[1]) * 1) |> as("dgCMatrix")
  x_mat_angfact_blunt <- ((x_mat_ang >= l_angles$critical_angles[2]) * 1) |> as("dgCMatrix")
  
  invisible(gc())
  message("test")
  gib <- list(
    x_sharp = x_mat_angfact_sharp,
    x_blunt = x_mat_angfact_blunt,
    l_angles = l_angles,
    data_info = list(
      samples = s_names,
      path_to_df_ang = path_to_write_angles
      # COMMENT: remove l_df_ang_md5... (which is the file name of the stored data.table)
      # name_df_ang_sharp = l_df_ang_md5[["sharp"]],
      # name_df_ang_blunt = l_df_ang_md5[["blunt"]]
    )
  )
  return(gib)
}