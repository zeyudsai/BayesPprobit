% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{sample_latent}
\alias{sample_latent}
\title{Sample Latent Variables}
\usage{
sample_latent(X, theta, y, p)
}
\arguments{
\item{X}{Design matrix}

\item{theta}{Current parameter values}

\item{y}{Binary response vector}

\item{p}{Shape parameter}
}
\value{
Vector of sampled latent variables
}
\description{
Samples latent variables from truncated p-GGD distribution
}
\keyword{internal}
