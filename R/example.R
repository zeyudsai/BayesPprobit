#=========================================
# BayesGProbit Package Basic Usage Example
#=========================================

library(BayesGProbit)

#------------------
# 1. Simulate Data
#------------------
set.seed(123)

# Generate sample dataset
n <- 5000  # number of observations
d <- 5     # number of features

# Create design matrix
X <- matrix(rnorm(n*d), n, d)
colnames(X) <- paste0("X", 1:d)

# True parameters
beta_true <- c(1.5, -0.8, 0.6, -0.4, 0.2)
p_true <- 2

# Generate response
eta <- X %*% beta_true
prob <- pgnorm(eta, alpha = p_scale(p_true), beta = p_true)
y <- rbinom(n, 1, prob)

#------------------------
# 2. Basic Model Fitting
#------------------------

# Fit model with default settings
fit_basic <- multi_chain(
  n_sim = 2000,    # number of iterations
  burn_in = 500,   # burn-in period
  X = X,
  y = y,
  initial_theta = rep(0, d),
  p = 2,          # initial p value
  n_chains = 3    # number of chains
)

fit_basic
# Print basic summary
summary(fit_basic)

#------------------------
# 3. Diagnostic Plots
#------------------------

# Create diagnostic plots
par(mfrow = c(2,2))

# Trace plot for p
plot(fit_basic$p_chains[[1]], type = "l",
     xlab = "Iteration", ylab = "p",
     main = "Trace Plot of p Parameter")
abline(h = p_true, col = "red", lty = 2)

# Posterior density of p
hist(fit_basic$p_chains[[1]], freq = FALSE,
     main = "Posterior Distribution of p",
     xlab = "p")
abline(v = p_true, col = "red", lty = 2)

# Compare true vs estimated parameters
plot(beta_true, fit_basic$posterior_beta,
     xlab = "True Values", ylab = "Estimated Values",
     main = "True vs Estimated Parameters")
abline(0, 1, col = "red", lty = 2)

# Reset plot layout
par(mfrow = c(1,1))

#------------------------
# 4. Using Coresets
#------------------------

# Generate larger dataset
n_large <- 50000
X_large <- matrix(rnorm(n_large*d), n_large, d)
y_large <- rbinom(n_large, 1, 0.5)

# Compare different coreset methods
methods <- c("leverage", "uniform", "oneshot")
coreset_size <- 1000

results <- list()
for(method in methods) {
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
    X = X_large[cs$indices,],
    y = y_large[cs$indices],
    initial_theta = rep(0, d),
    p = 2,
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

#------------------------
# 6. Visualization
#------------------------

# Plot posterior distributions of p for different methods
colors <- c("black", "red", "blue", "green")
plot(density(fit_basic$p_chains[[1]]),
     main = "Posterior Distribution of p",
     xlab = "p", col = colors[1])
lines(density(results$leverage$p_chains[[1]]), col = colors[2])
lines(density(results$uniform$p_chains[[1]]), col = colors[3])
lines(density(results$oneshot$p_chains[[1]]), col = colors[4])
legend("topright",
       legend = c("Full Data", "Leverage", "Uniform", "Oneshot"),
       col = colors, lty = 1)

#------------------------
# 7. Timing Comparison
#------------------------

# Compare computation times
times <- data.frame(
  Method = c("Full Data", "Leverage", "Uniform", "Oneshot"),
  Runtime = c(
    fit_basic$runtime[3],
    results$leverage$runtime[3],
    results$uniform$runtime[3],
    results$oneshot$runtime[3]
  )
)

print(times)

#------------------------
# 8. Prediction Example
#------------------------

# Generate test data
n_test <- 1000
X_test <- matrix(rnorm(n_test*d), n_test, d)

# Compute predicted probabilities
pred_prob <- pgnorm(
  X_test %*% fit_basic$posterior_beta,
  alpha = p_scale(fit_basic$posterior_p),
  beta = fit_basic$posterior_p
)

# Plot predicted probabilities
hist(pred_prob,
     main = "Distribution of Predicted Probabilities",
     xlab = "Probability", freq = FALSE)

#---------------------------------
# 9. Save Results for Future Use
#---------------------------------

# Create results directory if it doesn't exist
if(!dir.exists("results")) dir.create("results")

# Save model fits
saveRDS(fit_basic, "results/fit_basic.rds")
saveRDS(results, "results/coreset_results.rds")

#---------------------------------
# Session Information
#---------------------------------
sessionInfo()
