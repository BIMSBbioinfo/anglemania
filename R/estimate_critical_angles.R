#' Estimate critical angles
#'
#' @description
#' A wrapper around functions that estimate_critical_angles
#' and their levels of significance from an approximated
#' distribution of angles.
#'
#' @details
#' First, the angle distribution is approximated by quantile
#' splitting provided by the user. Then the gaussian is fit by
#' quantile transformation. Goodness of fit is tested by
#' Kolmogorov-Smirnov test. Measure of additional variance is
#' estimated from the data by integrating the difference between
#' the fitted gaussian and theoretical distribution. Critical
#' angles are extracted from the fitted gaussian.
#'
#' @importFrom magrittr %>%
#' @param x_mat_ang data.table. A long data.table with three columns: ########
#'   x - name of a gene in row, y - name of a gene in column,
#'   angle - an angle between x and y.
#' @param s_dims integer. Number of samples in the input gene
#'   expression matrix.
#' @param extrema double. Fraction of the angles
#'   to be cut from both sides of an approximated angle
#'   distribution.
#' @return list. Records on parameters of the gaussian,
#' goodness of fit, and critical angles.
#' @export estimate_critical_angles
estimate_critical_angles <- function(x_mat_ang, #nolint
                                     s_dims,
                                     extrema) {
  quantile_split <- min(extrema)
  ##
  out <- approximate_angles(x_mat_ang, quantile_split) %>%
    fit_gaussian(.) %>%
    test_gof_Ks(.) %>%
    estimate_MOAV(., s_dims) %>%
    get_analytical_extremes(., extrema)
  ##
  return(out)
}