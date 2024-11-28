# utils.R

#' @title Utility Functions for p-Generalized Probit Regression
#' @description Helper functions for working with p-generalized Gaussian distribution (p-GGD)
#' @name BayesPprobit-utils
#' @return A collection of utility functions:
#' \describe{
#'   \item{p_scale}{Returns the scale parameter value for p-GGD}
#'   \item{truncated_gnorm}{Returns samples from truncated p-GGD}
#' }
NULL

#' @title Scale Parameter for p-GGD
#' @description Computes scale parameter for p-generalized Gaussian distribution
#' @param p Shape parameter
#' @return A numeric value representing the scale parameter p^(1/p)
#' @export
p_scale <- function(p) {
  if (!is.numeric(p) || length(p) != 1 || p <= 0) {
    stop("`p` must be a single positive numeric value.")
  }
  p^(1 / p)
}

#' @rdname BayesPprobit-utils
#' @param n Number of samples to generate
#' @param range Vector of length 2 specifying the truncation range
#' @param mu Location parameter
#' @param alpha Scale parameter
#' @param beta Shape parameter
#' @return A numeric vector of samples from the truncated distribution
#' @keywords internal
truncated_gnorm <- function(n, range, mu, alpha, beta) {
  F_a <- gnorm::pgnorm(min(range), mu = mu, alpha = alpha, beta = beta)
  F_b <- gnorm::pgnorm(max(range), mu = mu, alpha = alpha, beta = beta)

  u <- runif(n, min = F_a, max = F_b)
  u[u < 1e-6] <- 1e-6
  u[u > 1 - 1e-6] <- 1 - 1e-6

  result <- gnorm::qgnorm(u, mu = mu, alpha = alpha, beta = beta)
  return(result)
}
