#' Fit a guassian distribution
#'
#' @description
#' Fits the gaussian distribution to the approximated angle distribution
#' using the sum of residual squares.
#'
#' @details
#' Fitting is done using non-linear minimization from the **stats** package.
#'
#' @importFrom tidyr tibble
#' @param l_approx list. An output from the **approxinate.angles** function.
#' @return list. List with multiple entries, one of which contains
#' an approxiamted distribution and details of the gaussian fit.
#' Others contains NAs to be filled by subsequent functions.
#' @export fit_gaussian
fit_gaussian <- function(l_approx) { #nolint
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