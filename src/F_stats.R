## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
estimate.critical.angles <- function(x_df_ang, #nolint
                                     s_dims, 
                                     extrema) {
  ##
  quantile_split <- min(extrema)
  ## --
  out <- approximate.angles(x_df_ang, quantile_split) %>%
    fit.gaussian(.) %>%
    test.gof.Ks(.) %>%
    estimate.MOAV(., s_dims) %>%
    get.analytical.extremes(., extrema)
  ##
  return(out)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
approximate.angles <- function(x_df_ang, #nolint
                               quantile_split) {
  ##
  if (quantile_split <= 0) {
    stop("percentile.split should be a positive number")
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
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
fit.gaussian <- function(l_approx) { #nolint
  ## function to compute sum of residual squares to a fit
  f <- function(q, x, prob) {
    res <- pnorm(x, q[1], q[2]) - prob
    sum(res * res)
  }
  ## fit the least squares and extract the coefficients
  coeff <- (fit <- nlm(f, c(90, 10),
                       l_approx$angles_dist$angle,
                       l_approx$angles_dist$prob))$estimate
  ## add estimated cumulative probabilities
  l_approx$angles_dist$prob_gauss <- pnorm(l_approx$angles_dist$angle, 
                                           mean = coeff[1],
                                           sd = coeff[2])
  ##
  l_approx$statistics$mu <- coeff[1]
  l_approx$statistics$sigma <- coeff[2]
  ##
  return(l_approx)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
test.gof.Ks <- function(l_approx, #nolint
                        alpha = 0.01,
                        mt_nsims = 6000,
                        seed = 42) {
  ##
  mu_est <- l_approx$statistics$mu
  sigma_est <- l_approx$statistics$sigma
  ref_angles <- l_approx$angles_dist$angle
  ref_probs <- l_approx$angles_dist$prob
  ##
  initial_d <- max(abs(ref_probs - l_approx$angles_dist$prob_gauss))
  background_d <- montecarlo.Dstat(ref_angles, ref_probs,
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
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## Monte-Carlo simulation of background D statistic distribution from
## estimated parameters of a Gaussian distribution
montecarlo.Dstat <- function(ref_data, #nolint
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
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
estimate.MOAV <- function(l_approx, #nolint
                          s_dims) {
  ##
  mu_est <- l_approx$statistics$mu
  sigma_est <- l_approx$statistics$sigma
  ## approximation of a PDF for a distribution of random angles on a sphere
  h <- function(theta, p) {
    (1 / sqrt(pi)) * (gamma(p / 2) / gamma((p - 1) / 2)) * sin(theta)^(p - 2)
  }
  ## simulate pdf
  df <- tidyr::tibble(
    angle_pi = seq(0, pi, length.out = 181),
    angle = 0:180,
    h_theta = h(seq(0, pi, length.out = 181), s_dims),
    gauss_est = dnorm(0:180, mean = mu_est, sd = sigma_est)
  ) %>%
    dplyr::mutate(h_theta_est = h_theta / sum(h_theta))
  ##
  moav <- df %>%
    dplyr::filter(gauss_est > h_theta_est) %>%
    dplyr::summarise(moav = sum(gauss_est - h_theta_est)) %>%
    dplyr::pull(moav)
  ##
  l_approx$angles_analytical <- df
  l_approx$statistics$moav <- moav
  l_approx$statistics$s_dims <- s_dims
  ##
  return(l_approx)
}
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
get.analytical.extremes <- function(l_approx, #nolint
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
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## merge angle statistics
extract.angle.stats <- function(l_processed) {
  purrr::map2_dfr(
    l_processed[["l_angles"]],
    names(l_processed[["l_angles"]]),
    ~ cbind(
      tibble(
        sample = .y,
        crit_sharp = .x$critical_angles[1],
        crit_blunt = .x$critical_angles[2]
      ),
      .x$statistics
    )
  )
}
