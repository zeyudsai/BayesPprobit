% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{multi_chain}
\alias{multi_chain}
\title{Multi-chain MCMC Sampling for Bayesian p-Generalized Probit Model}
\usage{
multi_chain(
  n_sim,
  burn_in,
  X,
  y,
  initial_theta,
  true_theta = NULL,
  initial_p = 2,
  mh_iter = 100,
  p_range = c(0.5, 5),
  step_size = 0.05,
  n_chains = 5
)
}
\arguments{
\item{n_sim}{Number of MCMC iterations}

\item{burn_in}{Number of burn-in iterations to discard}

\item{X}{Design matrix}

\item{y}{Binary response vector (0/1)}

\item{initial_theta}{Initial parameter values}

\item{true_theta}{True parameter values (optional, for simulation studies)}

\item{initial_p}{Initial value for the p parameter}

\item{mh_iter}{Number of Metropolis-Hastings iterations}

\item{p_range}{Range for p parameter sampling, vector of (min, max)}

\item{step_size}{Step size for MH proposals}

\item{n_chains}{Number of parallel chains}
}
\value{
A list containing:
\describe{
\item{beta_chains}{List of MCMC chains for beta parameters}
\item{p_chains}{List of chains for p parameter}
\item{posterior_beta}{Posterior mean estimates of beta}
\item{posterior_p}{Posterior mean of p}
\item{runtime}{Computation time}
\item{gelman_diag}{Gelman-Rubin convergence diagnostics}
\item{true_theta}{True parameter values (if provided)}
}
}
\description{
Performs multiple chain MCMC sampling for Bayesian p-generalized
probit regression using Metropolis-Hastings within Gibbs sampling.
}
