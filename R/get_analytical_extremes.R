#' Test goodness of fit
#'
#' @description
#' Test the goodness of fit of the gaussian with the Kolmogorov-Smirnov
#' test.
#'
#' @details
#' Goodness of fit is determined by comparing the analztical D against
#' the simulation background distribution of D. One-sided p-value is
#' computed agaisnt a provided alpha.
#'
#' @importFrom purrr map2_dfr
#' @param l_approx list. An output from the **approximate_angles** function.
#' @param extrema double. Fraction of the angles
#'   to be cut from both sides of an approximated angle
#'   distribution.
# @param alpha double. Significance level for the goodness of fit.
# @param mt_nsims integer. Number of simulated Ds.
# @param seed integer. Random seed. Important for reproducibility
#'   of the background simulation.
#' @return list. List with multiple entries, one of which contains
#' an approxiamted distribution and details of the gaussian fit.
#' Others contains NAs to be filled by subsequent functions.
#' @export get_analytical_extremes
get_analytical_extremes <- function(l_approx, #nolint
                                    extrema) {
  ##
  ref_angles <- l_approx$angles_dist$angle
  gauss_quants <- l_approx$angles_dist$prob_gauss
  ##
  crit_vals <- purrr::map_dbl(extrema,
                              ~ ref_angles[which.min(abs(gauss_quants - .x))]
  )
  ##
  l_approx$critical_angles <- crit_vals
  ##
  return(l_approx)
}