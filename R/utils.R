#' @title Scale Parameter for p-Generalized Gaussian Distribution
#'
#' @description
#' Computes the scale parameter for the p-generalized Gaussian distribution
#' to maintain unit variance.
#'
#' @param p Numeric. Shape parameter.
#'
#' @return Numeric. Scale parameter value.
#'
#' @details
#' The scale parameter is computed as p^(1/p) to ensure unit variance
#' of the p-GGD distribution.
#'
#' @keywords internal
p_scale <- function(p) {
  p^(1/p)
}

#' @title Log Likelihood for p-Generalized Probit Model
#'
#' @description
#' Computes the log likelihood for p-generalized probit regression.
#'
#' @param X Design matrix.
#' @param beta Regression coefficients.
#' @param p Shape parameter.
#' @param y Binary response vector.
#'
#' @return Numeric. Log likelihood value.
#'
#' @details
#' Computes:
#' sum(y * log(pi) + (1-y) * log(1-pi))
#' where pi = Phi_p(X*beta)
#'
#' Includes numerical safeguards for extreme probabilities.
#'
#' @keywords internal

#' @title p-Generalized Gaussian Distribution Functions
#' @description Density, distribution, quantile and random generation functions
#' for the p-generalized Gaussian distribution
#' @name pggd
NULL

#' @rdname pggd
#' @export
dgnorm <- function(x, mu = 0, alpha = 1, beta = 2) {
  # PDF of p-generalized Gaussian distribution
  z <- (abs(x - mu)/alpha)^beta
  beta/(2 * alpha * gamma(1/beta)) * exp(-z)
}

#' @rdname pggd
#' @export
pgnorm <- function(x, mu = 0, alpha = 1, beta = 2) {
  # CDF of p-generalized Gaussian distribution

  # For numerical stability
  if(length(x) > 1) {
    result <- sapply(x, function(xi) {
      pgnorm(xi, mu, alpha, beta)
    })
    return(result)
  }

  # Handle edge cases
  if(is.infinite(x)) {
    return(ifelse(x > 0, 1, 0))
  }

  # Calculate standardized value
  z <- (x - mu)/alpha

  # Split calculation for positive and negative z
  if(z >= 0) {
    0.5 + 0.5 * pgamma((z^beta)/beta, shape = 1/beta)
  } else {
    0.5 - 0.5 * pgamma(((-z)^beta)/beta, shape = 1/beta)
  }
}

#' @rdname pggd
#' @export
qgnorm <- function(p, mu = 0, alpha = 1, beta = 2) {
  # Quantile function of p-generalized Gaussian distribution

  # For vectors
  if(length(p) > 1) {
    result <- sapply(p, function(pi) {
      qgnorm(pi, mu, alpha, beta)
    })
    return(result)
  }

  # Handle edge cases
  if(p <= 0) return(-Inf)
  if(p >= 1) return(Inf)

  # For p < 0.5, use symmetry
  if(p < 0.5) {
    return(mu - alpha * (beta * qgamma(1 - 2*p, shape = 1/beta))^(1/beta))
  }

  # For p >= 0.5
  mu + alpha * (beta * qgamma(2*p - 1, shape = 1/beta))^(1/beta)
}

#' @rdname pggd
#' @export
rgnorm <- function(n, mu = 0, alpha = 1, beta = 2) {
  # Random generation for p-generalized Gaussian distribution
  u <- runif(n)
  qgnorm(u, mu, alpha, beta)
}

#' @title Scale Parameter for p-GGD
#' @description Computes scale parameter for p-generalized Gaussian distribution
#' @param p Shape parameter
#' @return Scale parameter value
#' @export
p_scale <- function(p) {
  (gamma(3/p)/gamma(1/p))^(1/2)
}
#p_scale <- function(p) {
#  p^(1/p)
#}

#' @title Log Likelihood for p-Probit Model
#' @description Computes log likelihood for p-generalized probit regression
#' @keywords internal
log_likelihood <- function(X, beta, p, y) {
  Xb <- X %*% beta
  pi <- pgnorm(Xb, alpha = p_scale(p), beta = p)

  # Avoid numerical issues
  pi[pi == 0] <- .Machine$double.xmin
  pi[pi == 1] <- 1 - .Machine$double.eps

  sum(y * log(pi) + (1-y) * log(1-pi))
}

#' @title Sample from Truncated p-GGD
#' @description Generates samples from truncated p-generalized Gaussian distribution
#' @keywords internal
rtrunc_pggd <- function(n, mu, p, lower = -Inf, upper = Inf) {
  # Implementation using inverse CDF method
  F.a <- pgnorm(lower, mu = mu, alpha = p_scale(p), beta = p)
  F.b <- pgnorm(upper, mu = mu, alpha = p_scale(p), beta = p)

  u <- runif(n, min = F.a, max = F.b)
  qgnorm(u, mu = mu, alpha = p_scale(p), beta = p)
}
