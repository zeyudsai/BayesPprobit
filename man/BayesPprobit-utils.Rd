% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{BayesPprobit-utils}
\alias{BayesPprobit-utils}
\alias{truncated_gnorm}
\title{Utility Functions for p-Generalized Probit Regression}
\usage{
truncated_gnorm(n, range, mu, alpha, beta)
}
\arguments{
\item{n}{Number of samples to generate}

\item{range}{Vector of length 2 specifying the truncation range}

\item{mu}{Location parameter}

\item{alpha}{Scale parameter}

\item{beta}{Shape parameter}
}
\value{
A collection of utility functions:
\describe{
\item{p_scale}{Returns the scale parameter value for p-GGD}
\item{truncated_gnorm}{Returns samples from truncated p-GGD}
}

A numeric vector of samples from the truncated distribution
}
\description{
Helper functions for working with p-generalized Gaussian distribution (p-GGD)
}
\keyword{internal}
