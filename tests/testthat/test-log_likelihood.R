# tests/testthat/test-log_likelihood.R

test_that("log_likelihood function works correctly", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  beta <- c(0.5, -0.5)
  p <- 2
  y <- c(1, 0, 1)

  pi <- pgnorm(X %*% beta, alpha = p_scale(p), beta = p)

  manual_ll <- sum(y * log(pi) + (1 - y) * log(1 - pi))

  function_ll <- log_likelihood(X, beta, p, y)

  expect_equal(function_ll, manual_ll, tolerance = 1e-6)

  expect_error(log_likelihood(X, beta, p = -1, y), "p must be a single positive numeric value.")
  expect_error(log_likelihood(X, beta, p = 2, y = c(1, 0)), "length(y) must match number of rows in X.")
})
