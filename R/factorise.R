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
#' @importFrom anglemana extract.angles melt.to.df
#' @importFrom anglemana estimate.critical.angles filter.angles write.angles
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
#' @export
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
  x_mat_ang <- extract.angles(x_mat, n_threads = n_threads)
  invisible(gc())
  message("Melting to data table...")
  x_df_ang <- melt.to.df(x_mat_ang)
  invisible(gc())
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
  l_flippity <- filter.angles(x_df_ang, l_angles)
  if (!is.na(path_to_write_angles) & is.character(path_to_write_angles)) {
    message("Writing tables with angles to compute integration metric...")
    l_df_ang_md5 <- write.angles(l_flippity, path_to_write_angles)
  }
  message("Factorising angle matrices...")
  x_mat_angfact_sharp <- ifelse(x_mat_ang <= l_angles$critical_angles[1], 1, 0)
  x_mat_angfact_sharp <- as(x_mat_angfact_sharp, "dgCMatrix")
  x_mat_angfact_blunt <- ifelse(x_mat_ang >= l_angles$critical_angles[2], 1, 0)
  x_mat_angfact_blunt <- as(x_mat_angfact_blunt, "dgCMatrix")
  invisible(gc())
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