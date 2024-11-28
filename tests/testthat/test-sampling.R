test_that("multi_chain works correctly", {
  set.seed(123)
  n <- 20
  d <- 3
  X <- matrix(rnorm(n * d), n, d)
  y <- rbinom(n, 1, 0.5)

  fit <- multi_chain(
    n_sim = 50,
    burn_in = 10,
    X = X,
    y = y,
    initial_theta = rep(0, d),
    initial_p = 2,
    mh_iter = 10,
    p_range = c(1, 3),
    step_size = 0.1,
    n_chains = 2
  )


  stopifnot(inherits(fit, "pgprobit"))
  stopifnot(all(c("beta_chains", "p_chains", "posterior_beta",
                  "posterior_p", "runtime", "gelman_diag",
                  "true_theta") %in% names(fit)))

  stopifnot(length(fit$beta_chains) == 2)
  stopifnot(length(fit$p_chains) == 2)
  stopifnot(length(fit$posterior_beta) == d)
  stopifnot(is.numeric(fit$posterior_p))
})

test_that("multi_chain handles errors", {
  n <- 20
  d <- 3
  X <- matrix(rnorm(n * d), n, d)
  y <- rbinom(n, 1, 0.5)

  expect_error(multi_chain(n_sim = 50, burn_in = 10, X = as.data.frame(X),
                           y = y, initial_theta = rep(0, d)))
  expect_error(multi_chain(n_sim = 50, burn_in = 10, X = X, y = 1:n,
                           initial_theta = rep(0, d)))
})
