## Test the prior constructors and the validation

test_that("prior_normal creates a valid prior object", {
  pr <- prior_normal(0, 2.5)

  expect_s3_class(pr, "BayesPprobit_prior")
  expect_equal(pr$family, "normal")
  expect_equal(pr$location, 0)
  expect_equal(pr$scale, 2.5)
  expect_equal(pr$df, 1)
})

test_that("prior_student_t creates a valid prior object", {
  pr <- prior_student_t(4, 0, 3)

  expect_s3_class(pr, "BayesPprobit_prior")
  expect_equal(pr$family, "student_t")
  expect_equal(pr$location, 0)
  expect_equal(pr$scale, 3)
  expect_equal(pr$df, 4)
})

test_that("prior_laplace creates a valid prior object", {
  pr <- prior_laplace(location = 0, scale = 1)

  expect_s3_class(pr, "BayesPprobit_prior")
  expect_equal(pr$family, "laplace")
  expect_equal(pr$location, 0)
  expect_equal(pr$scale, 1)
  expect_equal(pr$df, 1)
})

test_that("invalid prior parameters are rejected", {
  expect_error(prior_normal(scale = -1))
  expect_error(prior_normal(scale = 0))
  expect_error(prior_student_t(df = 0))
  expect_error(prior_student_t(df = -1))
  expect_error(prior_student_t(scale = -1))
  expect_error(prior_student_t(scale = 0))
  expect_error(prior_laplace(scale = -1))
  expect_error(prior_laplace(scale = 0))
})

## test stan_p_probit

test_that("stan_p_probit rejects non-binary outcomes", {
  X <- matrix(rnorm(20), ncol = 2)
  y <- c(rep(0, 9), 2)
  expect_error(
    stan_p_probit(X, y)
  )
})

test_that("stan_p_probit rejects invalid p bounds", {
  X <- matrix(rnorm(20), ncol = 2)
  y <- rbinom(10, 1, 0.5)

  expect_error(stan_p_probit(X, y, p_bounds = c(0, 6)))
  expect_error(stan_p_probit(X, y, p_bounds = c(3, 2)))
  expect_error(stan_p_probit(X, y, p_bounds = 2))
  expect_error(stan_p_probit(X, y, p_bounds = c(NA, 6)))
  expect_error(stan_p_probit(X, y, p_bounds = c(1, Inf)))
})

test_that("stan_p_probit rejects incompatible X and y", {
  X <- matrix(rnorm(20), ncol = 2)
  y <- rbinom(9, 1, 0.5)

  expect_error(stan_p_probit(X, y))
})

test_that("stan_p_probit rejects objects that are not BayesPprobit priors", {
  X <- matrix(rnorm(20), ncol = 2)
  y <- rbinom(10, 1, 0.5)

  expect_error(
    stan_p_probit(
      X,
      y,
      prior_beta = list(family = "normal", location = 0, scale = 2.5, df = 1)
    )
  )
})

test_that("stan_p_probit rejects unsupported prior families", {
  X <- matrix(rnorm(20), ncol = 2)
  y <- rbinom(10, 1, 0.5)

  bad_prior <- structure(
    list(
      family = "exponential",
      location = 0,
      scale = 1,
      df = 1
    ),
    class = "BayesPprobit_prior"
  )

  expect_error(stan_p_probit(
    X,
    y,
    prior_beta = bad_prior
    ))
})

test_that("stan_p_probit rejects missing values", {
  X <- matrix(rnorm(20), ncol = 2)
  y <- rbinom(10, 1, 0.5)

  X[1, 1] <- NA
  expect_error(stan_p_probit(X, y))

  X[1, 1] <- 0
  y[1] <- NA
  expect_error(stan_p_probit(X, y))
})

test_that("stan_p_probit runs with a small dataset", {
  set.seed(123)

  N <- 50
  X <- matrix(rnorm(N * 2), ncol = 2)
  eta <- -0.5 + X %*% c(0.5, -1)
  y <- rbinom(N, 1, pnorm(eta))

  fit <- stan_p_probit(
    X = X,
    y = y,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 200,
    seed = 123,
    refresh = 0
  )
  expect_true(inherits(fit, "CmdStanMCMC"))
})
