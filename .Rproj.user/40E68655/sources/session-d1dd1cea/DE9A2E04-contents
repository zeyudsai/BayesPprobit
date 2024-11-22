# tests/testthat/test-distributions.R

test_that("dgnorm returns correct density values", {
  x <- 0
  mu <- 0
  alpha <- 1
  beta <- 2

  expect_equal(dgnorm(x, mu, alpha, beta), dnorm(x, mean = mu, sd = alpha), tolerance = 1e-4)

  x <- c(-1, 0, 1)
  expected <- dnorm(x, mean = mu, sd = alpha)
  expect_equal(dgnorm(x, mu, alpha, beta), expected, tolerance = 1e-4)
})

test_that("pgnorm returns correct cumulative probabilities", {
  x <- 0
  mu <- 0
  alpha <- 1
  beta <- 2

  expect_equal(pgnorm(x, mu, alpha, beta), pnorm(x, mean = mu, sd = alpha), tolerance = 1e-4)

  x <- c(-1, 0, 1)
  expected <- pnorm(x, mean = mu, sd = alpha)
  expect_equal(pgnorm(x, mu, alpha, beta), expected, tolerance = 1e-4)
})

test_that("qgnorm returns correct quantiles", {
  p <- 0.5
  mu <- 0
  alpha <- 1
  beta <- 2

  expect_equal(qgnorm(p, mu, alpha, beta), qnorm(p, mean = mu, sd = alpha), tolerance = 1e-4)

  p <- c(0.025, 0.5, 0.975)
  expected <- qnorm(p, mean = mu, sd = alpha)
  expect_equal(qgnorm(p, mu, alpha, beta), expected, tolerance = 1e-4)
})

test_that("rgnorm generates random samples with correct properties", {
  set.seed(123)
  samples <- rgnorm(1000, mu = 0, alpha = 1, beta = 2)

  expect_length(samples, 1000)
  expect_true(abs(mean(samples)) < 0.1)
  expect_true(abs(sd(samples) - 1) < 0.1)
})
