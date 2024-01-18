#' Simulate background D
#'
#' @description
#' Monte-Carlo simulation of background D statistic distribution from
#' estimated parameters of a Gaussian distribution
#'
#' @details
#' The function samples **nsamp** times from a gaussian distribution
#' with mean and standard deviation provided by the user. Then, the
#' samples are tested against the provided reference and resulting
#' D-statistics are recorded in the output vector.
#'
#' @import tidyr
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @param ref_data vector double. Vector with angles
#'   corresponding to the probabilities in **ref_quants**.
#' @param ref_quants vector double. Vector with
#'   probabilities corresponding to the angles in **ref_data**.
#' @param mu double. Mean of the fitted gaussian.
#' @param sigma double. Standard deviation of the
#'   fitted gaussian.
#' @param nsamp integer. Number of the simulated distributions.
#' @param seed integer. Random seed. Important for reproducibility
#'   of the background simulation.
#' @return vector double. Vector of Kolmogorov-Smirnov D-statistics.
#' @export montecarlo_Dstat
montecarlo_Dstat <- function(ref_data, #nolint
                             ref_quants,
                             mu,
                             sigma,
                             nsamp = 2000,
                             seed = 42) {
  set.seed(seed)
  n <- length(ref_data)
  y <- rnorm(n * nsamp, mean = mu, sd = sigma)
  d <- tidyr::tibble(
    sim = rep(1:nsamp, each = n),
    y = y
  ) %>%
    tidyr::nest(data = y) %>%
    dplyr::mutate(
      ecdf = purrr::map(data, function(df) ecdf(df$y)),
      D = purrr::map(ecdf,
        function(f_ecdf) {
          test_quants <- f_ecdf(ref_data)
          max(abs(ref_quants - test_quants))
        }
      )
    ) %>%
    tidyr::unnest(D) %>%
    dplyr::pull(D)
  ##
  return(d)
}