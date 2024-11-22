# tests/testthat/test-methods.R

test_that("print.pgprobit works correctly", {
  # 创建一个模拟 pgprobit 对象
  pgprobit_obj <- list(
    beta_chains = list(matrix(rnorm(100), nrow = 10)),
    p_chains = list(runif(100)),
    posterior_beta = rnorm(10),
    posterior_p = runif(1),
    runtime = 1.23,
    gelman_diag = list(psrf = matrix(1.0, nrow = 10, ncol = 2)),
    true_theta = rnorm(10)
  )
  class(pgprobit_obj) <- "pgprobit"

  # 捕获打印输出
  expect_output(print(pgprobit_obj), "Bayesian p-Generalized Probit Regression")
  expect_output(print(pgprobit_obj), "Posterior mean of p:")
  expect_output(print(pgprobit_obj), "Posterior means of coefficients:")
  expect_output(print(pgprobit_obj), "Gelman-Rubin convergence diagnostics:")
})

test_that("plot.pgprobit runs without errors", {
  # 创建一个模拟 pgprobit 对象
  pgprobit_obj <- list(
    beta_chains = list(matrix(rnorm(100), nrow = 10)),
    p_chains = list(runif(100)),
    posterior_beta = rnorm(10),
    posterior_p = runif(1),
    runtime = 1.23,
    gelman_diag = list(psrf = matrix(1.0, nrow = 10, ncol = 2)),
    true_theta = rnorm(10)
  )
  class(pgprobit_obj) <- "pgprobit"

  # 测试绘图函数运行不出错
  expect_silent(plot(pgprobit_obj))
})
