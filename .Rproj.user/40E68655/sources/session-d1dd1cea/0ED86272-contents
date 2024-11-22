# tests/testthat/test-distributions.R
test_that("dgnorm returns correct density values", {
  x <- 0
  mu <- 0
  beta <- 2
  alpha <- 1 * sqrt(2)  # Adjust alpha to match sd = 1 in dnorm

  actual <- dgnorm(x, mu = mu, alpha = alpha, beta = beta)
  expected <- dnorm(x, mean = mu, sd = 1)

  expect_equal(actual, expected)
})
test_that("pgnorm returns correct cumulative probabilities", {
  x <- c(-1, 0, 1)
  mu <- 0
  beta <- 2
  alpha <- 1 * sqrt(2)

  actual <- pgnorm(x, mu = mu, alpha = alpha, beta = beta)
  expected <- pnorm(x, mean = mu, sd = 1)

  expect_equal(actual, expected)
})

test_that("qgnorm returns correct quantiles", {
  p <- c(0.025, 0.5, 0.975)
  mu <- 0
  beta <- 2
  alpha <- 1 * sqrt(2)

  actual <- qgnorm(p, mu = mu, alpha = alpha, beta = beta)
  expected <- qnorm(p, mean = mu, sd = 1)

  expect_equal(actual, expected)
})

test_that("rgnorm generates random samples with correct properties", {
  set.seed(123)
  n <- 10000
  mu <- 0
  beta <- 2
  alpha <- 1 * sqrt(2)

  samples <- rgnorm(n, mu = mu, alpha = alpha, beta = beta)
  sample_mean <- mean(samples)
  sample_sd <- sd(samples)

  expect_true(abs(sample_mean - mu) < 0.05)
  expect_true(abs(sample_sd - 1) < 0.05)
})

