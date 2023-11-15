#' Estiamte additional variance
#'
#' @description
#' Estimates the measure of additional variance (MOAV) in reference
#' to the expected varience.
#'
#' @details
#' The function integrates non-overlapping area between the
#' analytical and theoritical distributions of angles between
#' vectors of gene expression.
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom tidyr tibble
#' @param l_approx list. List with angular statistics after the
#'   fitting of a gaussian.
#' @param s_dims integer. Vector with
#'   probabilities corresponding to the angles in **ref_data**.
#' @return vector double. Vector of Kolmogorov-Smirnov D-statistics.
estimate.MOAV <- function(l_approx, #nolint
                          s_dims) {
  mu_est <- l_approx$statistics$mu
  sigma_est <- l_approx$statistics$sigma
  ## approximation of a PDF for a distribution of random angles on a sphere
  h <- function(theta, p) {
    (1 / sqrt(pi)) * (gamma(p / 2) / gamma((p - 1) / 2)) * sin(theta)^(p - 2)
  }
  df <- tidyr::tibble(
    angle_pi = seq(0, pi, length.out = 181),
    angle = 0:180,
    h_theta = h(seq(0, pi, length.out = 181), s_dims),
    gauss_est = dnorm(0:180, mean = mu_est, sd = sigma_est)
  ) %>%
    dplyr::mutate(h_theta_est = h_theta / sum(h_theta))
  moav <- df %>%
    dplyr::filter(gauss_est > h_theta_est) %>%
    dplyr::summarise(moav = sum(gauss_est - h_theta_est)) %>%
    dplyr::pull(moav)
  l_approx$angles_analytical <- df
  l_approx$statistics$moav <- moav
  l_approx$statistics$s_dims <- s_dims
  return(l_approx)
}