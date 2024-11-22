# sampling.R

#' @title Multi-chain MCMC Sampling for Bayesian p-Generalized Probit Model
#' @description Performs multiple chain MCMC sampling for Bayesian p-generalized
#'   probit regression using Metropolis-Hastings within Gibbs sampling.
#' @param n_sim Number of MCMC iterations
#' @param burn_in Number of burn-in iterations to discard
#' @param X Design matrix
#' @param y Binary response vector (0/1)
#' @param initial_theta Initial parameter values
#' @param true_theta True parameter values (optional, for simulation studies)
#' @param p Initial value for the p parameter
#' @param mh_iter Number of Metropolis-Hastings iterations
#' @param p_range Range for p parameter sampling, vector of (min, max)
#' @param step_size Step size for MH proposals
#' @param n_chains Number of parallel chains
#' @return A list containing:
#'   \item{beta_chains}{List of MCMC chains}
#'   \item{p_chains}{List of chains for p parameter}
#'   \item{posterior_beta}{Posterior mean estimates}
#'   \item{posterior_p}{Posterior mean of p}
#'   \item{runtime}{Computation time}
#'   \item{gelman_diag}{Gelman-Rubin convergence diagnostics}
#'   \item{true_theta}{True parameter values (if provided)}
#' @export
multi_chain <- function(n_sim,
                        burn_in,
                        X,
                        y,
                        initial_theta,
                        true_theta = NULL,
                        p = 2,
                        mh_iter = 100,
                        p_range = c(0.1, 5),
                        step_size = 0.01,
                        n_chains = 5) {

  # Input validation
  if(!methods::is(X, "matrix")) stop("X must be a matrix")
  if(!is.numeric(y) || !all(y %in% c(0,1)))
    stop("y must be binary (0/1)")
  if(length(y) != nrow(X))
    stop("Length of y must match rows of X")
  if(length(initial_theta) != ncol(X))
    stop("Length of initial_theta must match cols of X")

  # Initialize storage
  beta_chains <- list()
  p_chains <- list()
  start_time <- proc.time()

  # Create progress bar
  pb <- utils::txtProgressBar(min = 0, max = n_chains, style = 3)

  # Run multiple chains
  for(i in 1:n_chains) {
    result <- gibbs_sampler(
      n_sim = n_sim,
      burn_in = burn_in,
      X = X,
      y = y,
      initial_theta = initial_theta,
      p = p,
      mh_iter = mh_iter,
      p_range = p_range,
      step_size = step_size
    )

    beta_chains[[i]] <- result$beta_chain
    p_chains[[i]] <- result$p_chain

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  # Compute runtime
  runtime <- proc.time() - start_time

  # Convert chains to mcmc objects for diagnostics
  mcmc_chains <- coda::mcmc.list(lapply(beta_chains, coda::mcmc))

  # Compute convergence diagnostics
  gelman_diag <- coda::gelman.diag(mcmc_chains)

  # Compute summary statistics
  post_beta <- matrix(0, nrow=n_chains, ncol=length(initial_theta))
  post_p <- numeric(n_chains)

  for(i in 1:n_chains) {
    post_beta[i,] <- colMeans(beta_chains[[i]][-(1:burn_in),])
    post_p[i] <- mean(p_chains[[i]][-(1:burn_in)])
  }

  # Prepare result
  result <- list(
    beta_chains = beta_chains,
    p_chains = p_chains,
    posterior_beta = colMeans(post_beta),
    posterior_p = mean(post_p),
    runtime = runtime,
    gelman_diag = gelman_diag,
    true_theta = true_theta
  )

  # Assign S3 class
  class(result) <- "pgprobit"

  # Return results
  return(result)
}


#' @title Sample Latent Variables
#' @keywords internal
sample_latent <- function(X, theta, y, p) {
  n <- length(y)
  mu <- as.vector(X %*% theta)
  alpha <- p_scale(p)

  # Initialize latent variables
  z <- numeric(n)

  # For y = 0, sample from (-∞, 0)
  idx_0 <- which(y == 0)
  if(length(idx_0) > 0) {
    z[idx_0] <- truncnorm::rtruncnorm(
      n = length(idx_0),
      mean = mu[idx_0],
      sd = alpha,
      a = -Inf,
      b = 0
    )
  }

  # For y = 1, sample from (0, ∞)
  idx_1 <- which(y == 1)
  if(length(idx_1) > 0) {
    z[idx_1] <- truncnorm::rtruncnorm(
      n = length(idx_1),
      mean = mu[idx_1],
      sd = alpha,
      a = 0,
      b = Inf
    )
  }

  return(z)
}

#' @title Update Parameters via Gibbs Step
#' @keywords internal
update_theta <- function(X, z, p) {
  # Get dimensions
  n <- nrow(X)
  d <- ncol(X)

  # Compute scale parameter
  alpha <- p_scale(p)

  # Solve linear system for mean
  XtX <- Matrix::crossprod(X)
  Xtz <- Matrix::crossprod(X, z)
  mu <- Matrix::solve(XtX, Xtz)

  # Compute covariance matrix
  Sigma <- alpha^2 * Matrix::solve(XtX)

  # Draw from multivariate normal
  as.vector(mvtnorm::rmvnorm(1, mean = mu, sigma = Sigma))
}

#' @title Update Shape Parameter via MH
#' @keywords internal
update_p <- function(X, z, theta, p, mh_iter, p_range, step_size) {
  p_current <- p

  # Compute current log-likelihood
  ll_current <- sum(stats::dnorm(
    z - as.vector(X %*% theta),
    mean = 0,
    sd = p_scale(p_current),
    log = TRUE
  ))

  for(i in 1:mh_iter) {
    # Propose new p value
    p_prop <- stats::runif(1,
                           min = max(p_range[1], p_current - step_size),
                           max = min(p_range[2], p_current + step_size))

    # Compute proposed log-likelihood
    ll_prop <- sum(stats::dnorm(
      z - as.vector(X %*% theta),
      mean = 0,
      sd = p_scale(p_prop),
      log = TRUE
    ))

    # Accept/reject
    log_ratio <- ll_prop - ll_current
    if(log(stats::runif(1)) < log_ratio) {
      p_current <- p_prop
      ll_current <- ll_prop
    }
  }

  return(p_current)
}

#' @title Gibbs Sampler Implementation
#' @keywords internal
gibbs_sampler <- function(n_sim, burn_in, X, y, initial_theta,
                          p, mh_iter, p_range, step_size) {
  # Get dimensions
  n <- nrow(X)
  d <- ncol(X)

  # Initialize
  theta <- initial_theta
  beta_chain <- matrix(0, nrow = n_sim, ncol = d)
  p_chain <- numeric(n_sim)

  # Store initial values
  beta_chain[1,] <- theta
  p_chain[1] <- p

  # Create progress bar
  pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)

  # Main loop
  for(t in 2:n_sim) {
    # Sample latent variables
    z <- sample_latent(X, theta, y, p)

    # Update theta
    theta <- update_theta(X, z, p)

    # Update p
    p <- update_p(X, z, theta, p, mh_iter, p_range, step_size)

    # Store results
    beta_chain[t,] <- theta
    p_chain[t] <- p

    # Update progress bar
    utils::setTxtProgressBar(pb, t)
  }
  close(pb)

  list(
    beta_chain = beta_chain,
    p_chain = p_chain
  )
}

#' @title Print Method for pgprobit Objects
#' @description Prints a summary of the Bayesian p-Generalized Probit Regression results.
#' @param x An object of class \code{pgprobit}.
#' @param ... Additional arguments passed to or from other methods.
#' @export
print.pgprobit <- function(x, ...) {
  cat("Bayesian p-Generalized Probit Regression\n\n")
  cat("Posterior mean of p:", round(x$posterior_p, 4), "\n")
  cat("Posterior means of coefficients:\n")
  print(round(x$posterior_beta, 4))
  cat("\nGelman-Rubin convergence diagnostics:\n")
  print(x$gelman_diag)
}

#' @title Plot Method for pgprobit Objects
#' @description Plots diagnostics for Bayesian p-Generalized Probit Regression results.
#' @param x An object of class \code{pgprobit}.
#' @param ... Additional arguments passed to or from other methods.
#' @export
plot.pgprobit <- function(x, ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  # Set up plotting layout
  graphics::par(mfrow = c(2, 2))

  # Trace plot of p
  graphics::plot(x$p_chains[[1]], type = "l",
                 xlab = "Iteration", ylab = "p",
                 main = "Trace Plot of p")

  # Density plot of p
  dens <- stats::density(x$p_chains[[1]])
  graphics::plot(dens, main = "Posterior Density of p",
                 xlab = "p")

  # Trace plots of some beta coefficients
  graphics::plot(x$beta_chains[[1]][,1], type = "l",
                 xlab = "Iteration", ylab = expression(beta[1]),
                 main = "Trace Plot of β1")

  # Gelman-Rubin statistics
  grstat <- x$gelman_diag$psrf[,1]
  graphics::barplot(grstat,
                    names.arg = paste("β", 1:length(grstat), sep=""),
                    main = "Gelman-Rubin Statistics",
                    ylab = "PSRF")
  graphics::abline(h = 1.1, col = "red", lty = 2)
}
