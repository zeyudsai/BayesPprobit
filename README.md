```markdown
# BayesPProbit

<!-- badges: start -->
[![R-CMD-check](https://github.com/username/BayesPProbit/workflows/R-CMD-check/badge.svg)](https://github.com/username/BayesPProbit/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/BayesPProbit)](https://CRAN.R-project.org/package=BayesPProbit)
<!-- badges: end -->

## Overview

**BayesPProbit** is an advanced R package that implements Bayesian p-generalized probit regression models, incorporating efficient data compression techniques through the use of coresets. This package is designed to handle large-scale binary classification problems by providing scalable and efficient Bayesian inference methods.

Key features of **BayesPProbit** include:

- **Flexible Modeling**: Supports the p-generalized probit model, allowing for robustness to outliers and heavy-tailed data by adjusting the shape parameter \( p \).
- **Efficient MCMC Sampling**: Implements Metropolis-Hastings within Gibbs sampling for posterior inference, optimized for performance.
- **Coreset Integration**: Incorporates coreset construction methods to reduce data size while preserving essential statistical properties, enabling scalable inference on massive datasets.
- **Diagnostic Tools**: Provides functions for convergence diagnostics, posterior summary statistics, and visualization of MCMC results.
- **Easy-to-Use Interface**: Designed with a user-friendly API that integrates seamlessly with existing R workflows.

## Installation

You can install the stable version of **BayesPProbit** from CRAN:

```r
install.packages("BayesPProbit")
```

Or install the development version from GitHub to access the latest features and updates:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install BayesPProbit from GitHub
devtools::install_github("username/BayesPProbit")
```

*Note*: Replace `"username"` with the actual GitHub repository owner if necessary.

## Usage

Below is a basic example demonstrating how to use **BayesPProbit** for Bayesian p-generalized probit regression:

```r
library(BayesPProbit)

# Simulate some data
set.seed(123)
n <- 1000  # Number of observations
d <- 5     # Number of predictors
X <- matrix(rnorm(n * d), n, d)
beta_true <- c(1.5, -0.8, 0.6, -0.4, 0.2)
p_true <- 2
eta <- X %*% beta_true
alpha <- p_scale(p_true)
prob <- gnorm::pgnorm(eta, mu = 0, alpha = alpha, beta = p_true)
y <- rbinom(n, 1, prob)

# Fit the Bayesian p-generalized probit model
fit <- multi_chain(
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

# Print a summary of the results
print(fit)

# Plot diagnostic plots
plot(fit)
```

For more detailed examples and advanced usage, please refer to the package vignette:

```r
vignette("BayesPProbit")
```

## Features

- **Bayesian Inference for Probit Models**: Implements a Bayesian approach to probit regression, allowing uncertainty quantification and incorporation of prior information.
- **p-Generalized Distribution**: Utilizes the p-generalized normal distribution to model the latent variables, offering flexibility in capturing data characteristics.
- **Metropolis-Hastings within Gibbs Sampling**: Employs efficient MCMC algorithms for sampling from complex posterior distributions.
- **Coreset Methods**: Supports various coreset construction techniques such as leverage score sampling, uniform sampling, and oneshot methods to handle large datasets.
- **Convergence Diagnostics**: Integrates with the `coda` package for convergence assessment using Gelman-Rubin diagnostics and other MCMC diagnostic tools.
- **Visualization**: Provides functions for trace plots, posterior density plots, and comparison of true vs. estimated parameters.

## Documentation

Comprehensive documentation is available for all functions within the package. Access the documentation using the standard R help system:

```r
# For a general overview
?BayesPProbit

# For specific functions
?multi_chain
?gibbs_sampler
```

## Support and Contributions

If you encounter any issues, have questions, or would like to contribute to the development of **BayesPProbit**, please visit our GitHub repository:

[GitHub - username/BayesPProbit](https://github.com/username/BayesPProbit)

Feel free to open issues or submit pull requests. Contributions are welcome!

## Citation

If you use **BayesPProbit** in your research or publications, please cite it as follows:

```r
citation("BayesPProbit")
```

This will provide the appropriate citation information, including authors and version number.

## License

**BayesPProbit** is licensed under the MIT License. See the [LICENSE](https://github.com/username/BayesPProbit/blob/master/LICENSE) file for more details.

## Acknowledgments

We would like to thank all contributors and users who have provided feedback and suggestions to improve the package.

---

*This README was generated to provide a comprehensive overview of the **BayesPProbit** package, its features, and how to get started. For more detailed information, please refer to the package documentation and vignettes.*
```

---

**Note**: Replace `"username"` with the actual GitHub username or organization that maintains the **BayesPProbit** package. Additionally, ensure that all links and references are updated to point to the correct locations in your repository.