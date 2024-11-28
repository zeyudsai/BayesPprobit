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
#' @param initial_p Initial value for the p parameter
#' @param mh_iter Number of Metropolis-Hastings iterations
#' @param p_range Range for p parameter sampling, vector of (min, max)
#' @param step_size Step size for MH proposals
#' @param n_chains Number of parallel chains
#' @return A list containing:
#'   \describe{
#'     \item{beta_chains}{List of MCMC chains for beta parameters}
#'     \item{p_chains}{List of chains for p parameter}
#'     \item{posterior_beta}{Posterior mean estimates of beta}
#'     \item{posterior_p}{Posterior mean of p}
#'     \item{runtime}{Computation time}
#'     \item{gelman_diag}{Gelman-Rubin convergence diagnostics}
#'     \item{true_theta}{True parameter values (if provided)}
#'   }
#' @export
#'
multi_chain <- function(n_sim,
                        burn_in,
                        X,
                        y,
                        initial_theta,
                        true_theta = NULL,
                        initial_p = 2,
                        mh_iter = 100,
                        p_range = c(0.5, 5),
                        step_size = 0.05,
                        n_chains = 5) {
  # Input validation
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!is.numeric(y) || !all(y %in% c(0, 1)))
    stop("y must be binary (0/1)")
  if (length(y) != nrow(X))
    stop("Length of y must match rows of X")
  if (length(initial_theta) != ncol(X))
    stop("Length of initial_theta must match cols of X")

  # Initialize storage
  beta_chains <- list()
  p_chains <- list()
  start_time <- proc.time()

  # Run multiple chains
  for (i in 1:n_chains) {
    result <- gibbs_sampler(
      n_sim = n_sim,
      burn_in = burn_in,
      X = X,
      y = y,
      initial_theta = initial_theta,
      initial_p = initial_p,
      mh_iter = mh_iter,
      p_range = p_range,
      step_size = step_size
    )

    beta_chains[[i]] <- result$beta_chain
    p_chains[[i]] <- result$p_chain
  }

  # Compute runtime
  runtime <- proc.time() - start_time

  # Convert chains to mcmc objects for diagnostics
  mcmc_chains <- coda::mcmc.list(lapply(beta_chains, function(chain)
    coda::mcmc(chain[-(1:burn_in), ])))

  # Compute convergence diagnostics
  gelman_diag <- coda::gelman.diag(mcmc_chains)

  # Compute summary statistics
  post_beta <- matrix(0, nrow = n_chains, ncol = length(initial_theta))
  post_p <- numeric(n_chains)

  for (i in 1:n_chains) {
    post_beta[i, ] <- colMeans(beta_chains[[i]][-(1:burn_in), ])
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

  return(result)
}

#' @title Sample Latent Variables
#' @description Samples latent variables from truncated p-GGD distribution
#' @param X Design matrix
#' @param theta Current parameter values
#' @param y Binary response vector
#' @param p Shape parameter
#' @return Vector of sampled latent variables
#' @keywords internal
sample_latent <- function(X, theta, y, p) {
  n <- length(y)
  mu <- as.vector(X %*% theta)

  # Initialize latent variables
  z <- numeric(n)

  # Compute alpha and beta parameters for gnorm functions
  alpha <- p_scale(p)
  beta <- p

  # For y = 1, sample from truncated distribution on (0, infinity)
  idx_1 <- which(y == 1)
  if (length(idx_1) > 0) {
    z[idx_1] <- truncated_gnorm(
      n = length(idx_1),
      range = c(0, Inf),
      mu = mu[idx_1],
      alpha = alpha,
      beta = beta
    )
  }

  # For y = 0, sample from truncated distribution on (-infinity, 0)
  idx_0 <- which(y == 0)
  if (length(idx_0) > 0) {
    z[idx_0] <- truncated_gnorm(
      n = length(idx_0),
      range = c(-Inf, 0),
      mu = mu[idx_0],
      alpha = alpha,
      beta = beta
    )
  }

  return(z)
}

#' @title Update Parameters via Gibbs Step
#' @description Updates regression coefficients using Gibbs sampling
#' @param X Design matrix
#' @param z Latent variables
#' @param p Shape parameter
#' @return Updated parameter vector
#' @keywords internal
update_theta <- function(X, z, p) {
  # Use sirt::lq_fit to estimate coefficients
  model <- sirt::lq_fit(y = z, X = X, pow = p,
                        est_pow = FALSE)
  beta_hat <- model$coefficients

  # Compute Sigma = tau(p) * (X^T X)^{-1}
  XtX <- crossprod(X)
  XtX_inv <- solve(XtX)
  tau_p <- p^(2 / p) * gamma(3 / p) / gamma(1 / p)
  Sigma <- tau_p * XtX_inv

  # Sample from multivariate normal distribution
  # theta_new <- mvtnorm::rmvnorm(1, mean = beta_hat, sigma = Sigma)
  theta_new <- numeric(ncol(X))
  for (i in 1:ncol(X)) {
    theta_new[i] <- gnorm::rgnorm(1, mu=beta_hat[i],
                       alpha =Sigma[i,i],
                       beta = p)
  }
  return(as.vector(theta_new))
}

#' @title Update Shape Parameter via Metropolis-Hastings
#' @description Updates p parameter using Metropolis-Hastings step
#' @param X Design matrix
#' @param z Latent variables
#' @param theta Current parameters
#' @param p Current p value
#' @param mh_iter Number of MH iterations
#' @param p_range Range for p values
#' @param step_size Step size for proposals
#' @return Updated p value
#' @keywords internal
#' @title Adaptive Step Size Metropolis-Hastings
#' @export
update_p <- function(X, z, theta, p, mh_iter, p_range, step_size) {
  Xb <- X%*%theta
  ehat <- z - Xb
  u <- runif(mh_iter)
  tol <- 1e-40
  p_current <- p
  accept_count <- 0
  target_rate <- 0.234

  step_history <- numeric(mh_iter)

  for (i in 1:mh_iter) {
    p_prop <- runif(1,
                    min = max(p_range[1], p_current - step_size),
                    max = min(p_range[2], p_current + step_size))

    R0 <- sum(log(gnorm::dgnorm(ehat,alpha = p_scale(p_prop),
                         beta = p_prop)+tol))
    R1 <- sum(log(dgnorm(ehat,alpha = p_scale(p_current),
                         beta = p_current)+tol))
    adj<- max(R0,R1)
    R  <- exp(R0-adj)/exp(R1-adj)


    if (u[i] < R) {
      p_current <- p_prop
      accept_count <- accept_count + 1
    }

    # if (i > 50) {
    #   current_rate <- accept_count / i
    #   if (current_rate < target_rate) {
    #     step_size <- step_size * 0.9
    #   } else {
    #     step_size <- step_size * 1.1
    #   }
    # }

    step_history[i] <- step_size
  }

  return(p_current)
}


#' @title Gibbs Sampler Implementation
#' @description Internal implementation of Gibbs sampler for p-generalized probit model
#' @param n_sim Number of iterations
#' @param burn_in Burn-in period
#' @param X Design matrix
#' @param y Response vector
#' @param initial_theta Initial parameters
#' @param initial_p Initial p value
#' @param mh_iter MH iterations
#' @param p_range Range for p
#' @param step_size Step size
#' @return List containing:
#'   \describe{
#'     \item{beta_chain}{MCMC chain for beta parameters}
#'     \item{p_chain}{Chain for p parameter}
#'   }
#' @keywords internal
gibbs_sampler <- function(n_sim, burn_in, X, y, initial_theta, initial_p,
                          mh_iter, p_range, step_size, fix_p = FALSE) {
  n <- nrow(X)
  d <- ncol(X)

  # Initialize
  theta <- initial_theta
  p <- initial_p
  beta_chain <- matrix(0, nrow = n_sim, ncol = d)
  p_chain <- numeric(n_sim)

  # Store initial values
  beta_chain[1, ] <- theta
  p_chain[1] <- p

  # Create progress bar
  pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)

  # Main loop
  for (t in 2:n_sim) {
    # Sample latent variables
    z <- sample_latent(X, theta, y, p)

    # Update p
    if(fix_p){
      p <- initial_p
    } else{
      p <- update_p(X, z, theta, p, mh_iter,
                    p_range, step_size)
    }



    # Update theta
    theta <- update_theta(X, z, p)

    # Store results
    beta_chain[t, ] <- theta
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
#' @param x An object of class \code{pgprobit}
#' @param ... Additional arguments passed to or from other methods
#' @return No return value, called for side effects of printing summary information.
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
#' @description Creates diagnostic plots for Bayesian p-Generalized Probit Regression results.
#' @param x An object of class \code{pgprobit}
#' @param ... Additional arguments passed to or from other methods
#' @return No return value, called for side effects of creating diagnostic plots.
#' @export
plot.pgprobit <- function(x, ...) {
  # Save old par settings
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
                 main = "Trace Plot of beta1")

  # Gelman-Rubin statistics
  grstat <- x$gelman_diag$psrf[,1]
  graphics::barplot(grstat,
                    names.arg = paste("beta", 1:length(grstat), sep=""),
                    main = "Gelman-Rubin Statistics",
                    ylab = "PSRF")
  graphics::abline(h = 1.1, col = "red", lty = 2)
}

