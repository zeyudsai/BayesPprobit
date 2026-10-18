# BayesPprobit

<!-- badges: start -->
[![R-CMD-check](https://github.com/zeyudsai/BayesPprobit/workflows/R-CMD-check/badge.svg)](https://github.com/zeyudsai/BayesPprobit/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/BayesPprobit)](https://CRAN.R-project.org/package=BayesPprobit)
<!-- badges: end -->

## Overview

**BayesPprobit** is an advanced R package that implements Bayesian p-generalized probit regression models, incorporating efficient data compression techniques through the use of coresets. This package is designed to handle large-scale binary classification problems by providing scalable and efficient Bayesian inference methods.

Key features of **BayesPprobit** include:

- **Flexible Modeling**: Supports the p-generalized probit model, allowing for robustness to outliers and heavy-tailed data by adjusting the shape parameter \( p \).
- **Flexible MCMC Sampling**: Provides Metropolis-Hastings within Gibbs sampling as well as Stan-based inference using the No-U-Turn Sampler (NUTS), an adaptive Hamiltonian Monte Carlo algorithm.
- **Coreset Integration**: Incorporates coreset construction methods to reduce data size while preserving essential statistical properties, enabling scalable inference on massive datasets.
- **Diagnostic Tools**: Provides functions for convergence diagnostics, posterior summary statistics, and visualization of MCMC results.
- **Easy-to-Use Interface**: Designed with a user-friendly API that integrates seamlessly with existing R workflows.

## Installation

You can install the stable version of **BayesPprobit** from CRAN:
```r
install.packages("BayesPprobit")
```

> **Note:** The CRAN release is currently pending final unit tests and edge case validation. The package will be submitted to CRAN upon completion. Please use the following installation.

Or install the development version from GitHub to access the latest features and updates:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install BayesPprobit from GitHub
devtools::install_github("zeyudsai/BayesPprobit")
```

### Stan backend

The Stan-based implementation requires CmdStanR and a working installation of CmdStan.

If CmdStan is not already installed, it can be installed using:

```r
cmdstanr::install_cmdstan()
```

## Usage

Below is a basic example demonstrating how to use **BayesPprobit** for Bayesian p-generalized probit regression:

```r
library(BayesPprobit)

# Simulate some data
set.seed(123)
n <- 1000   # Number of observations
d <- 5      # Number of predictors
X <- matrix(rnorm(n * d), n, d)
beta_true <- c(1.5, -0.8, 0.6, -0.4, 0.2)
intercept_true <- 0
p_true <- 2
eta <- intercept_true + X %*% beta_true
scale_p <- p_scale(p_true)
prob <- gnorm::pgnorm(eta, mu = 0, alpha = scale_p, beta = p_true)
y <- rbinom(n, 1, prob)

# Fit the Bayesian p-generalized probit model using the Gibbs sampler

fit_gibbs <- multi_chain(
  n_sim = 5000,
  burn_in = 1000,
  X = X,
  y = y,
  initial_theta = rep(0, d),
  initial_p = 2,
  mh_iter = 100,
  p_range = c(1, 3),
  step_size = 0.05,
  n_chains = 3
)

# Fit the Bayesian p-generalized probit model using Stan/NUTS

fit_stan <- stan_p_probit(
  X = X,
  y = y,
  prior_intercept = prior_normal(0, 5),
  prior_beta = prior_student_t(4, 0, 3),
  p_bounds = c(1, 6),
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  seed = 123,
  max_treedepth = 12
)

# Summarize results
print(fit_gibbs)
fit_stan$summary()

# Stan diagnostics
fit_stan$diagnostic_summary()

# Plot Gibbs sampler output
plot(fit_gibbs)
```
For more detailed examples and advanced usage, please refer to the package vignettes:

```r
vignette("BayesPprobit")
vignette("stan-p-probit")
```

## Features

- **Bayesian Inference for Probit Models**: Implements a Bayesian approach to probit regression, allowing uncertainty quantification and incorporation of prior information.
- **p-Generalized Distribution**: Utilizes the p-generalized normal distribution to model the latent variables, offering flexibility in capturing data characteristics.
- **Metropolis-Hastings within Gibbs Sampling**: Employs efficient MCMC algorithms for sampling from complex posterior distributions.
- **Stan-Based Inference**: Supports posterior sampling with the No-U-Turn Sampler (NUTS), an adaptive Hamiltonian Monte Carlo algorithm, via CmdStanR.
- **Coreset Methods**: Supports various coreset construction techniques such as leverage score sampling, uniform sampling, and oneshot methods to handle large datasets.
- **Convergence Diagnostics**: Supports convergence assessment for the Gibbs sampler through `coda`, while Stan fits provide the diagnostics available through CmdStanR, including \(\hat R\), effective sample sizes, divergences, and HMC diagnostics.
- **Visualization**: Provides visualization tools for MCMC output.

## Documentation

Comprehensive documentation is available for all functions within the package. Access the documentation using the standard R help system:

```r
# For a general overview
?BayesPprobit

# For specific functions
?multi_chain
?gibbs_sampler
?stan_p_probit
?prior_normal
?prior_student_t
?prior_laplace
```

## Support and Contributions

If you encounter any issues, have questions, or would like to contribute to the development of **BayesPprobit**, please visit our GitHub repository:

[GitHub - zeyudsai/BayesPprobit](https://github.com/zeyudsai/BayesPprobit)

Feel free to open issues or submit pull requests. Contributions are welcome!

## Citation

If you use **BayesPprobit** in your research or publications, please cite it as follows:

```r
citation("BayesPprobit")
```

You can also visit [CITATION](https://github.com/zeyudsai/BayesPprobit/blob/main/inst/CITATION) This will provide the appropriate citation information, including authors and version number.

## License

**BayesPprobit** is licensed under the MIT License. See the [LICENSE](https://github.com/zeyudsai/BayesPprobit/blob/main/LICENSE) file for more details.

## Acknowledgments

We would like to thank all contributors and users who have provided feedback and suggestions to improve the package.

---

*This README was generated to provide a comprehensive overview of the **BayesPprobit** package, its features, and how to get started. For more detailed information, please refer to the package documentation and vignettes.*
```

---
