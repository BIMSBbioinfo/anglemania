#' Estimate critical angles
#'
#' @description
#' A wrapper around functions that estimate critical angles
#' and their levels of signficance from an approximated
#' distribution of angles.
#'
#' @details
#' First, the angle distribution is approximated by quantile
#' splitting providid by the user. Then the gaussian is fit by
#' quantile transformation. Goodness of fit is tested by
#' Kolmogorov-Smirnov test. Measure of additional variance is
#' estiated fron the data by integrating the difference between
#' the fitted gaussian and theoretical distribution. Critical
#' angles are extracted from the fitted gaussian.
#'
#' @importFrom magrittr %>%
#' @importFrom anglemana approximate.angles fit.gaussian test.gof.Ks
#' @importFrom anglemana estimate.MOAV get.analytical.extremes
#' @param x_df_ang data.table. A long data.table with three columns:
#'   x - name of a gene in row, y - name of a gene in column,
#'   angle - an angle between x and y.
#' @param s_dims integer. Number of samples in the input gene
#'   expression matrix.
#' @param extrema double. Fraction of the angles
#'   to be cut from both sides of an approximated angle
#'   distribution.
#' @return list. Records on parameters of the gaussian,
#' goodness of fit, and critical angles.
estimate.critical.angles <- function(x_df_ang, #nolint
                                     s_dims,
                                     extrema) {
  quantile_split <- min(extrema)
  ##
  out <- approximate.angles(x_df_ang, quantile_split) %>%
    fit.gaussian(.) %>%
    test.gof.Ks(.) %>%
    estimate.MOAV(., s_dims) %>%
    get.analytical.extremes(., extrema)
  ##
  return(out)
}