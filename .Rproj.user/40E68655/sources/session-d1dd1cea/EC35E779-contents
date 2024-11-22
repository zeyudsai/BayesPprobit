# BayesPProbit/R/utils.R

# Load necessary packages
# Ensure these packages are listed in the DESCRIPTION file under Imports
library(pgnorm)
library(mvtnorm)
library(Matrix)
library(truncnorm)
library(abind)

#' Scale parameter for p-generalized Gaussian distribution
#'
#' @param p Shape parameter
#' @return Scale parameter
#' @export
p_scale <- function(p){
  if (!is.numeric(p) || length(p) != 1 || p <= 0) {
    stop("`p` must be a single positive numeric value.")
  }
  # Adjusting the scale to ensure p=2 corresponds to standard normal
  result <- ifelse(p == 2, 1, (gamma(3/p) / gamma(1/p))^(1/2))
  return(result)
}

#' Generate samples from truncated p-generalized Gaussian distribution
#'
#' @param n Number of samples
#' @param range Vector of two values indicating the truncation range
#' @param mu Mean parameter
#' @param alpha Scale parameter
#' @param beta Shape parameter
#' @return Vector of samples
#' @export
truncated_gnorm <- function(n, range, mu, alpha, beta) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("`n` must be a single positive integer.")
  }
  if (!is.numeric(range) || length(range) != 2) {
    stop("`range` must be a numeric vector of length 2.")
  }
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("`mu` must be a single numeric value.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("`alpha` must be a single positive numeric value.")
  }
  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0) {
    stop("`beta` must be a single positive numeric value.")
  }

  F_a <- pgnorm::pgnorm(min(range), alpha = alpha, beta = beta, mu = mu)
  F_b <- pgnorm::pgnorm(max(range), alpha = alpha, beta = beta, mu = mu)

  u <- runif(n, min = F_a, max = F_b)
  u[u < 1e-6] <- 1e-6
  u[u > 1 - 1e-6] <- 1 - 1e-6

  result <- qgnorm_custom(u, mu = mu, alpha = alpha, beta = beta)
  return(result)
}

#' Custom quantile function for p-generalized Gaussian distribution
#'
#' @param u Probability values (between 0 and 1)
#' @param mu Mean parameter
#' @param alpha Scale parameter
#' @param beta Shape parameter
#' @param lower Lower bound for quantile search
#' @param upper Upper bound for quantile search
#' @param tol Tolerance for root finding
#' @param maxiter Maximum iterations for root finding
#' @return Quantiles corresponding to probabilities u
#' @export
qgnorm_custom <- function(u, mu, alpha, beta, lower = -Inf, upper = Inf, tol = 1e-8, maxiter = 1000) {
  if (!is.numeric(u) || any(u < 0) || any(u > 1)) {
    stop("All elements of `u` must be between 0 and 1.")
  }
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("`mu` must be a single numeric value.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("`alpha` must be a single positive numeric value.")
  }
  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0) {
    stop("`beta` must be a single positive numeric value.")
  }

  # Define a helper function to find the root for a single u
  find_quantile <- function(prob) {
    # Initial guesses
    if (prob < 0.5) {
      lower_bound <- mu - 10 * alpha
      upper_bound <- mu
    } else {
      lower_bound <- mu
      upper_bound <- mu + 10 * alpha
    }

    # Define the objective function
    obj_fun <- function(x) {
      pgnorm::pgnorm(x, alpha = alpha, beta = beta, mu = mu) - prob
    }

    # Use tryCatch to handle cases where uniroot fails
    quantile <- tryCatch({
      uniroot(obj_fun, lower = lower_bound, upper = upper_bound, tol = tol, maxiter = maxiter)$root
    }, error = function(e) {
      warning(sprintf("Failed to find quantile for probability %.4f. Returning NA.", prob))
      return(NA)
    })

    return(quantile)
  }

  # Vectorize the quantile finding
  quantiles <- sapply(u, find_quantile)
  return(quantiles)
}

#' Update latent variable z in probit model
#'
#' @param X Design matrix
#' @param theta Current parameter estimates
#' @param y Response vector
#' @param Lp Shape parameter
#' @return Updated latent variable z
#' @export
latent_z <- function(X, theta, y, Lp){
  N <- length(y)
  N1 <- sum(y)
  N0 <- N - N1
  z <- rep(0, N)

  mu_z <- X %*% theta
  z[y == 0] <- truncated_gnorm(N0, mu = mu_z[y == 0], range = c(-Inf, 0), alpha = p_scale(Lp), beta = Lp)
  z[y == 1] <- truncated_gnorm(N1, mu = mu_z[y == 1], range = c(0, Inf), alpha = p_scale(Lp), beta = Lp)

  return(z)
}

#' Proposal Function for Metropolis Algorithm
#'
#' This function generates a proposed Lp value within the specified range.
#'
#' @param Lp Current Lp value.
#' @param range Vector of two values indicating the proposal range.
#' @param L Step size for the proposal.
#'
#' @return Proposed Lp value.
#' @export
proposal <- function(Lp, range, L) {
  proposed <- runif(1, min = max(min(range), Lp - L), max = min(max(range), Lp + L))
  return(proposed)
}

#' Metropolis Algorithm for Updating Lp
#'
#' This function updates the shape parameter Lp using the Metropolis algorithm.
#'
#' @param X Design matrix.
#' @param M_iter Number of Metropolis iterations.
#' @param theta Current theta values.
#' @param Lp Current Lp value.
#' @param range Range for proposing new Lp values.
#' @param L Step size for the proposal.
#' @param z Latent variable vector.
#' @param tol Tolerance to avoid log(0).
#'
#' @return Updated Lp value.
#' @export
Metropolis <- function(X, M_iter, theta, Lp, range, L, z, tol = 1e-40) {
  Xb <- X %*% theta
  ehat <- z - Xb

  for (i in 1:M_iter) {
    # Propose a new Lp value
    proposed <- proposal(Lp, range, L)

    # Calculate log-likelihood for current and proposed Lp
    R0 <- sum(log(pgnorm::dgnorm(ehat, alpha = p_scale(proposed), beta = proposed) + tol))
    R1 <- sum(log(pgnorm::dgnorm(ehat, alpha = p_scale(Lp), beta = Lp) + tol))

    # Calculate acceptance ratio
    adj <- max(R0, R1)
    R <- exp(R0 - adj) / exp(R1 - adj)

    # Accept or reject the proposed Lp
    if (runif(1) < R) {
      Lp <- proposed
    }
  }

  return(Lp)
}

#' Metropolis-Lp Gibbs Sampler
#'
#' This function performs Gibbs sampling with Metropolis updates for the shape parameter Lp.
#'
#' @param N_sim Number of MCMC iterations.
#' @param burn_in Number of burn-in iterations.
#' @param X Design matrix.
#' @param y Response vector.
#' @param initial_theta Initial theta values.
#' @param true_theta True theta values (for comparison).
#' @param Lp Initial Lp value.
#' @param M_iter Number of Metropolis iterations for updating Lp.
#' @param range Range for proposing new Lp values.
#' @param step Step size for Metropolis proposal.
#'
#' @return A list containing:
#' \describe{
#'   \item{chain}{Matrix of theta samples.}
#'   \item{modeling.time}{Total computation time.}
#'   \item{true_theta}{True theta values.}
#'   \item{post.mean}{Posterior mean of theta after burn-in.}
#'   \item{Lp}{Vector of Lp samples.}
#'   \item{Lp.mean}{Posterior mean of Lp after burn-in.}
#' }
#' @export
Metropolis_Lp_gibbssampler <- function(N_sim, burn_in, X, y, initial_theta,
                                       true_theta, Lp, M_iter, range, step) {
  D <- ncol(X)
  N <- nrow(X)

  theta <- initial_theta
  z <- rep(0, N)
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  Lp_chain <- numeric(N_sim)

  theta_chain[1, ] <- theta
  Lp_chain[1] <- Lp

  # Start timing
  initialtime <- proc.time()

  for (t in 2:N_sim) {
    # Update latent variable z
    z <- latent_z(X, theta, y, Lp)

    # Fit the p-generalized Gaussian model
    model <- lq_fit(y = z, X = X, pow = Lp, est_pow = FALSE)
    M <- model$coefficients

    # Update theta using pgnorm's dgnorm (assuming it's a placeholder for actual sampling)
    # Here, it's assumed that you have a method to sample theta; adjust accordingly
    # For illustration, we'll use the posterior mean
    theta <- M

    # Store theta
    theta_chain[t, ] <- theta

    # Update Lp using Metropolis
    Lp <- Metropolis(
      X = X,
      M_iter = M_iter,
      theta = theta,
      Lp = Lp,
      range = range,
      L = step,
      z = z
    )

    # Store Lp
    Lp_chain[t] <- Lp

    # Progress output
    if (t %% floor(N_sim / 100) == 0) {
      cat(sprintf("Progress %d%%\n", floor(100 * t / N_sim)))
    }
  }

  # Calculate posterior means
  post.mean <- colMeans(theta_chain[-(1:burn_in), ])
  p.mean <- mean(Lp_chain[-(1:burn_in)])

  # End timing
  modeling.time <- proc.time() - initialtime

  # Compile results
  chainlist <- list(
    "chain" = theta_chain,
    "modeling.time" = modeling.time,
    "true_theta" = true_theta,
    "post.mean" = post.mean,
    "Lp" = Lp_chain,
    "Lp.mean" = p.mean
  )

  return(chainlist)
}

#' Compute Sketch for Coreset Construction
#'
#' @param X Design matrix
#' @param Lp Shape parameter
#' @param sketch_size Size of the sketch
#' @return Sketched matrix X_b
#' @export
Compute_sketch <- function(X, Lp, sketch_size){
  N <- nrow(X)
  D <- ncol(X)
  B_i <- sample(1:sketch_size, N, replace = TRUE)
  sigma_i <- sample(c(-1,1), N, replace = TRUE)
  X_b <- matrix(0, nrow = sketch_size, ncol = D)
  lambda_i <- rexp(N, 1)
  if(Lp != 2) sigma_i <- sigma_i / (lambda_i^(1/Lp))

  for (i in 1:N) {
    X_b[B_i[i], ] <- X_b[B_i[i], ] + sigma_i[i] * X[i, ]
  }
  return(X_b)
}

#' Uniform Coreset Sampling
#'
#' @param sketch_size Size of the sketch
#' @param coreset_size Desired coreset size
#' @param X Design matrix
#' @param Lp Shape parameter
#' @return Vector of sampled indices
#' @export
Uniform_coreset <- function(sketch_size, coreset_size, X, Lp){
  N <- nrow(X)
  D <- ncol(X)
  scores <- runif(N)
  sample_index <- sample(1:N, size = coreset_size, replace = FALSE, prob = scores)
  return(sample_index)
}

#' Leverage Score for Lp Norm
#'
#' @param X Design matrix
#' @param Lp Shape parameter
#' @return Vector of leverage scores
#' @export
lp_leverage <- function(X, Lp){
  qr_X <- qr(X)
  Q <- qr.Q(qr_X)
  leverage_scores <- apply(Q, 1, function(row) sum(abs(row)^Lp)^(1/Lp))
  return(leverage_scores)
}

#' Lp Coreset Sampling
#'
#' @param sketch_size Size of the sketch
#' @param coreset_size Desired coreset size
#' @param X Design matrix
#' @param Lp Shape parameter
#' @return Vector of sampled indices
#' @export
Lp_coreset <- function(sketch_size, coreset_size, X, Lp){
  N <- nrow(X)
  D <- ncol(X)
  scores <- lp_leverage(X, Lp = Lp)
  sample_index <- sample(1:N, size = coreset_size, replace = FALSE, prob = scores)
  return(sample_index)
}

#' Sensitivity Score Calculation
#'
#' @param sketch_size Size of the sketch
#' @param coreset_size Desired coreset size
#' @param X Design matrix
#' @param Lp Shape parameter
#' @return Vector of sensitivity scores
#' @export
Sensitivity_score <- function(sketch_size, coreset_size, X, Lp){
  N <- nrow(X)
  D <- ncol(X)
  X_b <- Compute_sketch(X, Lp = Lp, sketch_size = sketch_size)
  R_matrix <- qr.R(qr(X_b))

  if (Lp == 2 && log(N) < D) {
    G <- diag(rnorm(D, 0, 1 / log(N)))
  } else {
    G <- diag(D)
  }

  XR <- X %*% (solve(R_matrix) %*% G)
  qi <- apply(XR, 1, function(row) sum(abs(row)^Lp))

  scores <- qi
  return(scores)
}

#' One-shot Coreset Sampling Based on Sensitivity Scores
#'
#' @param sketch_size Size of the sketch
#' @param coreset_size Desired coreset size
#' @param X Design matrix
#' @param Lp_max Maximum shape parameter
#' @return Vector of combined sensitivity scores
#' @export
one_shot_coreset <- function(sketch_size, coreset_size, X, Lp_max){
  N <- nrow(X)
  delta <- 1 / log(N)
  score <- rep(0, N)
  l <- round(log(Lp_max) / log(1 + delta))

  for (i in 0:l) {
    p <- (1 + delta)^i
    score_i <- Sensitivity_score(sketch_size = sketch_size, coreset_size = coreset_size, X = X, Lp = p)
    score <- score + score_i
  }
  score <- score + 1 / N
  return(score)
}

#' One-shot Coreset Sampling with Multiple Shape Parameters
#'
#' @param sketch_size Size of the sketch
#' @param coreset_size Desired coreset size
#' @param X Design matrix
#' @param Lp_max Maximum shape parameter
#' @param times Number of shape parameter iterations
#' @return List containing combined sensitivity scores and modeling time
#' @export
one_shot_coreset_times <- function(sketch_size, coreset_size, X, Lp_max, times){
  initialtime <- proc.time()
  N <- nrow(X)
  delta <- 1 / log(N)
  score <- rep(0, N)

  for (i in 0:times) {
    p <- (1 + delta)^i
    score_i <- Sensitivity_score(sketch_size = sketch_size, coreset_size = coreset_size, X = X, Lp = p)
    score <- score + score_i
  }
  score <- score + 1 / N
  modeling.time <- proc.time() - initialtime

  coreset_list <- list("score" = score, "modeling.time" = modeling.time)
  return(coreset_list)
}

#' Gaussian Kernel Function
#'
#' @param x_1 Vector or matrix
#' @param x_2 Vector or matrix
#' @return Gaussian kernel value
#' @export
gaussian_kernel <- function(x_1, x_2){
  x_1 <- as.matrix(x_1)
  x_2 <- as.matrix(x_2)
  result <- exp(-0.5 * sum((x_1 - x_2)^2))
  return(result)
}

#' Polynomial Kernel Function
#'
#' @param x_1 Vector or matrix
#' @param x_2 Vector or matrix
#' @return Polynomial kernel value
#' @export
polynomial_kernel <- function(x_1, x_2){
  result <- (1 + sum(x_1 * x_2))^2
  return(result)
}

#' Maximum Mean Discrepancy (MMD) Calculation
#'
#' @param sample_x Vector or matrix of samples from distribution X
#' @param sample_y Vector or matrix of samples from distribution Y
#' @return MMD value
#' @export
mmd <- function(sample_x, sample_y){
  m <- length(sample_x)
  n <- length(sample_y)

  # Compute sum of Gaussian kernels for sample_x
  sum_x <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      sum_x <- sum_x + gaussian_kernel(sample_x[i], sample_x[j])
    }
  }

  # Compute sum of Gaussian kernels for sample_y
  sum_y <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      sum_y <- sum_y + gaussian_kernel(sample_y[i], sample_y[j])
    }
  }

  # Compute sum of Gaussian kernels between sample_x and sample_y
  sum_x_y <- 0
  for (i in 1:m) {
    for (j in 1:n) {
      sum_x_y <- sum_x_y + gaussian_kernel(sample_x[i], sample_y[j])
    }
  }

  # Calculate MMD
  result <- sqrt((1 / m^2) * sum_x + (1 / n^2) * sum_y - (2 / (m * n)) * sum_x_y)
  return(result)
}

#' Generate y from p-generalized probit model
#'
#' @param p Shape parameter
#' @param X Design matrix
#' @param Xb Linear predictor (X %*% beta)
#' @param N Number of observations
#' @return Vector of binary responses y
#' @export
generate_p_y <- function(p, Xb, N){
  pi <- pgnorm::pgnorm(Xb, alpha = p_scale(p), beta = p)
  y.p <- rbinom(N, 1, pi)
  return(y.p)
}
