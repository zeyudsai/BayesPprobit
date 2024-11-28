test_that("p_scale function works correctly", {
  diff1 <- abs(p_scale(2) - sqrt(2))
  stopifnot(diff1 < 1e-6)

  diff2 <- abs(p_scale(1) - 1)
  stopifnot(diff2 < 1e-6)

  diff3 <- abs(p_scale(3) - 3^(1/3))
  stopifnot(diff3 < 1e-6)


  expect_error(p_scale(-1))
  expect_error(p_scale(0))
  expect_error(p_scale(c(1,2)))
})

test_that("truncated_gnorm works correctly", {
  set.seed(123)
  n <- 100
  mu <- 0
  alpha <- 1
  beta <- 2

  samples_pos <- truncated_gnorm(n, c(0, Inf), mu, alpha, beta)
  stopifnot(min(samples_pos) >= 0)
  stopifnot(is.numeric(samples_pos))
  stopifnot(length(samples_pos) == n)

  samples_neg <- truncated_gnorm(n, c(-Inf, 0), mu, alpha, beta)
  stopifnot(max(samples_neg) <= 0)
  stopifnot(is.numeric(samples_neg))
  stopifnot(length(samples_neg) == n)
})
