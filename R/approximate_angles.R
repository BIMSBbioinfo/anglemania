#' Approximate angle distribution
#'
#' @description
#' Split the angle distribution by a vector of probabilities.
#'
#' @details
#' Vector of probabilities for a quantile function is constructed
#' from a minimal splitting unit provided by the user.
#'
#' @importFrom tidyr tibble
#' @param x_df_ang data.table. A long data.table with three columns:
#'   x - name of a gene in row, y - name of a gene in column,
#'   angle - an angle between x and y.
#' @param quantile_split double. minimal splitting unit to
#'   construct a vector of probabilities from 0 to 1.
#' @return list. List with multiple entries, one of which contains
#' an approxiamted distribution. Others contains NAs to be filled by
#' subsequent functions.
approximate_angles <- function(x_df_ang, #nolint
                               quantile_split) {
  ##
  if (quantile_split <= 0) {
    stop("quantile split should be a positive number")
  }
  angles_dist <- tidyr::tibble(
    angle = unname(
      quantile(x_df_ang$angle,
        seq(0, 1, by = quantile_split)
      )
    ),
    prob = seq(0, 1, by = quantile_split)
  )
  ##
  sts <- tidyr::tibble(
    modality = NA,
    mu = NA,
    sigma = NA,
    estimated_d = NA,
    critical_d = NA,
    alpha = NA,
    p_value = NA,
    s_dims = NA,
    moav = NA
  )
  ##
  l_approx <- list(critical_angles = NA,
                   angles_dist = angles_dist,
                   angles_anylitical = NA,
                   statistics = sts)
  ##
  return(l_approx)
}