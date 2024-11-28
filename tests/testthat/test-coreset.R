test_that("compute_coreset works correctly", {
  set.seed(123)
  n <- 100
  d <- 3
  X <- matrix(rnorm(n * d), n, d)
  y <- rbinom(n, 1, 0.5)

  cs <- compute_coreset(X, y, coreset_size = 20)


  stopifnot(is.list(cs))
  stopifnot(all(c("indices", "weights", "scores") %in% names(cs)))

  stopifnot(length(cs$indices) == 20)
  stopifnot(length(cs$weights) == 20)
  stopifnot(length(cs$scores) == n)

  stopifnot(is.numeric(cs$indices))
  stopifnot(is.numeric(cs$weights))
  stopifnot(is.numeric(cs$scores))
  stopifnot(all(cs$indices >= 1 & cs$indices <= n))
  stopifnot(all(cs$weights > 0))
  stopifnot(all(cs$scores >= 0))
})
