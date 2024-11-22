# tests/testthat/test-sampling.R

test_that("multi_chain function returns correct structure", {
  set.seed(123)
  X <- matrix(rnorm(100), nrow = 10, ncol = 10)
  y <- sample(0:1, 10, replace = TRUE)
  initial_theta <- rep(0, 10)

  result <- multi_chain(
    n_sim = 50,
    burn_in = 10,
    X = X,
    y = y,
    initial_theta = initial_theta,
    initial_p = 2,  # Correct argument
    mh_iter = 10,
    p_range = c(0.1, 5),
    step_size = 0.01,
    n_chains = 2
  )


  expect_type(result, "list")
  expect_length(result, 7)
  expect_named(result, c("beta_chains", "p_chains", "posterior_beta", "posterior_p", "runtime", "gelman_diag", "true_theta"))

  expect_s3_class(result, "pgprobit")

  expect_length(result$beta_chains, 2)
  expect_length(result$p_chains, 2)
  expect_true(all(sapply(result$beta_chains, is.matrix)))
  expect_true(all(sapply(result$p_chains, is.numeric)))

  expect_length(result$posterior_beta, 10)
  expect_type(result$posterior_p, "double")
})

