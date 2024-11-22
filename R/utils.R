# utils.R

#' @title p-Generalized Gaussian Distribution Functions
#' @description Density, distribution, quantile and random generation functions
#' for the p-generalized Gaussian distribution
#' @name pgnorm_functions
NULL

#' @title Scale Parameter for p-GGD
#' @description Computes scale parameter for p-generalized Gaussian distribution
#' @param p Shape parameter
#' @return Scale parameter value
#' @export
p_scale <- function(p) {
  if (!is.numeric(p) || length(p) != 1 || p <= 0) {
    stop("`p` must be a single positive numeric value.")
  }
  p^(1 / p)
}

#' @rdname pgnorm_functions
#' @param x Numeric vector of quantiles
#' @param mu Location parameter
#' @param p Shape parameter
#' @param scale Scale parameter (optional)
#' @param log Logical; if TRUE, return log density
#' @export
d_pgnorm <- function(x, mu = 0, p = 2, scale = NULL, log = FALSE) {
  if (p <= 0) {
    stop("Parameter p must be positive.")
  }

  if (is.null(scale)) {
    scale <- p_scale(p)
  }
  beta <- p
  alpha <- scale

  # Check for numerical issues with gamma function
  gamma_val <- gamma(1 / beta)
  if (!is.finite(gamma_val) || gamma_val == 0) {
    if (log) {
      return(rep(-Inf, length(x)))
    } else {
      return(rep(0, length(x)))
    }
  }

  coef <- beta / (2 * alpha * gamma_val)
  exponent <- - (abs(x - mu) / alpha)^beta

  if (log) {
    return(log(coef) + exponent)
  } else {
    return(coef * exp(exponent))
  }
}

#' @rdname pgnorm_functions
#' @param q Numeric vector of quantiles
#' @param lower.tail logical; if TRUE, probabilities are smaller than x
#' @export
p_pgnorm <- function(q, mu = 0, p = 2, scale = NULL, lower.tail = TRUE) {
  if (is.null(scale)) {
    scale <- p_scale(p)
  }

  z <- (q - mu) / scale

  cdf_values <- ifelse(
    z >= 0,
    0.5 + 0.5 * pgamma((z^p) / p, shape = 1 / p),
    0.5 - 0.5 * pgamma(((-z)^p) / p, shape = 1 / p)
  )

  if (!lower.tail) {
    cdf_values <- 1 - cdf_values
  }

  return(cdf_values)
}

#' @rdname pgnorm_functions
#' @param n Number of observations
#' @export
r_pgnorm <- function(n, mu = 0, p = 2, scale = NULL) {
  if (is.null(scale)) {
    scale <- p_scale(p)
  }
  u <- runif(n)
  q_pgnorm(u, mu = mu, p = p, scale = scale)
}

#' @rdname pgnorm_functions
#' @param p_vals Numeric vector of probabilities
#' @export
q_pgnorm <- function(p_vals, mu = 0, p = 2, scale = NULL) {
  if (is.null(scale)) {
    scale <- p_scale(p)
  }

  n <- length(p_vals)
  if(length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if(length(scale) == 1) {
    scale <- rep(scale, n)
  }

  sapply(seq_len(n), function(i) {
    p_val <- p_vals[i]
    mu_i <- mu[i]
    scale_i <- scale[i]

    lower_bound <- mu_i - 10 * scale_i
    upper_bound <- mu_i + 10 * scale_i

    # Ensure that lower_bound < upper_bound
    if (lower_bound >= upper_bound) {
      lower_bound <- mu_i - abs(mu_i) - 10
      upper_bound <- mu_i + abs(mu_i) + 10
    }

    tryCatch({
      uniroot(function(x) p_pgnorm(x, mu = mu_i, p = p, scale = scale_i) - p_val,
              interval = c(lower_bound, upper_bound))$root
    }, error = function(e) {
      warning(paste("Failed to find root for p =", p_val, "mu =", mu_i))
      NA
    })
  })
}

#' @rdname pgnorm_functions
#' @keywords internal
rtrunc_pgnorm <- function(n, mu = 0, p = 2, scale = NULL, lower = -Inf, upper = Inf) {
  if (is.null(scale)) {
    scale <- p_scale(p)
  }

  if(length(mu) == 1) mu <- rep(mu, n)
  if(length(lower) == 1) lower <- rep(lower, n)
  if(length(upper) == 1) upper <- rep(upper, n)

  F_lower <- p_pgnorm(lower, mu, p, scale)
  F_upper <- p_pgnorm(upper, mu, p, scale)

  u <- runif(n, F_lower, F_upper)

  q_pgnorm(u, mu, p, scale)
}

#' @rdname pgnorm_functions
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
