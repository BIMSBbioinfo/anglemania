#' Test goodness of fit
#'
#' @description
#' Kolmogorov-Smirnov test for the goodness of fit.
#'
#' @details
#' Goodness of fit is determined by comparing the analztical D against
#' the simulation background distribution of D. One-sided p-value is
#' computed against a provided alpha.
#'
##' @importFrom anglemania montecarlo_Dstat
#' @param l_approx list. An output from the **approxinate.angles** function.
#' @param alpha double. Significance level for the goodness of fit.
#' @param mt_nsims integer. The number of simulated distributions.
#' @param seed integer. Random seed. Important for reproducibility
#'   of the background simulation.
#' @return list. List with multiple entries, one of which contains
#' an approxiamted distribution and details of the gaussian fit.
#' Others contains NAs to be filled by subsequent functions.
#' @export test_gof_Ks
test_gof_Ks <- function(l_approx, #nolint
                        alpha = 0.01,
                        mt_nsims = 6000,
                        seed = 42) {
  mu_est <- l_approx$statistics$mu
  sigma_est <- l_approx$statistics$sigma
  ref_angles <- l_approx$angles_dist$angle
  ref_probs <- l_approx$angles_dist$prob
  ##
  initial_d <- max(abs(ref_probs - l_approx$angles_dist$prob_gauss))
  background_d <- montecarlo_Dstat(ref_angles, ref_probs,
                                   mu = mu_est,
                                   sigma = sigma_est,
                                   nsamp = mt_nsims,
                                   seed = seed)
  critical_d <- quantile(background_d, 1 - alpha)
  ##
  l_approx$statistics$estimated_d <- initial_d
  l_approx$statistics$critical_d <- critical_d
  l_approx$statistics$alpha <- alpha
  l_approx$statistics$p_value <- (
    sum(background_d[background_d > initial_d]) / sum(background_d)
  )
  ##
  return(l_approx)
}