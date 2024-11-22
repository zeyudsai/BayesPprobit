#' @title Example Usage of BayesPprobit Package
#' @description Demonstrates basic usage of the package through examples
#' @export
run_example <- function() {
  #------------------
  # 1. Simulate Data
  #------------------
  set.seed(123)

  # Generate sample dataset
  n <- 1000  # number of observations
  d <- 5     # number of features

  # Create design matrix
  X <- matrix(stats::rnorm(n * d), n, d)
  colnames(X) <- paste0("X", 1:d)

  # True parameters
  beta_true <- c(1.5, -0.8, 0.6, -0.4, 0.2)
  p_true <- 2

  # Generate response
  eta <- X %*% beta_true  # linear predictor
  alpha <- p_scale(p_true)  # scale parameter
  prob <- gnorm::pgnorm(eta, mu = 0,  alpha = alpha, beta = p_true)
  y <- stats::rbinom(n, 1, prob)

  #------------------------
  # 2. Basic Model Fitting
  #------------------------

  # Fit model with default settings
  fit_basic <- multi_chain(
    n_sim = 10000,    # number of iterations
    burn_in = 500,   # burn-in period
    X = X,
    y = y,
    initial_theta = rep(0, d),
    initial_p = 2,         # initial p value
    n_chains = 3    # number of chains
  )

  # Print basic summary
  print(fit_basic)

  #------------------------
  # 3. Diagnostic Plots
  #------------------------

  # Create diagnostic plots
  graphics::par(mfrow = c(2,2))

  # Trace plot for p
  graphics::plot(fit_basic$p_chains[[1]], type = "l",
                 xlab = "Iteration", ylab = "p",
                 main = "Trace Plot of p Parameter")
  graphics::abline(h = p_true, col = "red", lty = 2)

  # Posterior density of p
  graphics::hist(fit_basic$p_chains[[1]], freq = FALSE,
                 main = "Posterior Distribution of p",
                 xlab = "p")
  graphics::abline(v = p_true, col = "red", lty = 2)

  # Compare true vs estimated parameters
  graphics::plot(beta_true, fit_basic$posterior_beta,
                 xlab = "True Values", ylab = "Estimated Values",
                 main = "True vs Estimated Parameters")
  graphics::abline(0, 1, col = "red", lty = 2)

  # Reset plot layout
  graphics::par(mfrow = c(1,1))

  #------------------------
  # 4. Using Coresets
  #------------------------

  # Generate larger dataset
  n_large <- 50000
  X_large <- matrix(stats::rnorm(n_large * d), n_large, d)
  eta_large <- X_large %*% beta_true
  prob_large <- p_pgnorm(eta_large, mu = 0, p = p_true, scale = alpha)
  y_large <- stats::rbinom(n_large, 1, prob_large)

  # Compare different coreset methods
  methods <- c("leverage", "uniform", "oneshot")
  coreset_size <- 1000

  results <- list()
  for (method in methods) {
    # Construct coreset
    cs <- compute_coreset(
      X = X_large,
      y = y_large,
      coreset_size = coreset_size,
      method = method
    )

    # Fit model on coreset
    fit_cs <- multi_chain(
      n_sim = 2000,
      burn_in = 500,
      X = X_large[cs$indices, ],
      y = y_large[cs$indices],
      initial_theta = rep(0, d),
      initial_p = 2,
      n_chains = 3
    )

    results[[method]] <- fit_cs
  }

  #------------------------
  # 5. Compare Results
  #------------------------

  # Compare posterior means across methods
  comparison <- data.frame(
    Parameter = c("p", paste0("beta", 1:d)),
    Full_Data = c(fit_basic$posterior_p, fit_basic$posterior_beta),
    Leverage = c(results$leverage$posterior_p, results$leverage$posterior_beta),
    Uniform = c(results$uniform$posterior_p, results$uniform$posterior_beta),
    Oneshot = c(results$oneshot$posterior_p, results$oneshot$posterior_beta)
  )

  print(comparison)

  # Return results invisibly
  invisible(list(
    basic_fit = fit_basic,
    coreset_results = results,
    comparison = comparison
  ))
}

#' @title Load Example Results
#' @description Loads pre-computed example results if available
#' @export
load_example_results <- function() {
  example_file <- system.file("examples", "example_results.rds",
                              package = "BayesPprobit")
  if (file.exists(example_file)) {
    readRDS(example_file)
  } else {
    stop("Example results not found. Run run_example() first.")
  }
}

# Run the example
#run_example()
