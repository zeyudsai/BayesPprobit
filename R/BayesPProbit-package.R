#' @title Bayesian p-Generalized Probit Regression
#'
#' @description
#' Implements Bayesian inference for p-generalized probit regression models with
#' efficient data compression using coresets. The package provides comprehensive
#' tools for:
#' \itemize{
#'   \item Fitting Bayesian p-generalized probit models via MCMC
#'   \item Efficient handling of large datasets through coreset construction
#'   \item Posterior inference and convergence diagnostics
#'   \item Visualization of MCMC results
#' }
#'
#' @details
#' The key features of the package include:
#' \itemize{
#'   \item Flexible p-generalized probit models that extend beyond standard probit regression
#'   \item Efficient MCMC sampling using Metropolis-Hastings within Gibbs
#'   \item Coreset construction methods for scalable computation with large datasets
#'   \item Multiple chain sampling with convergence diagnostics
#'   \item Visualization tools for posterior analysis
#' }
#'
#' The main functions are:
#' \describe{
#'   \item{multi_chain}{Run MCMC with multiple chains}
#'   \item{compute_coreset}{Construct coresets for data compression}
#'   \item{plot.pgprobit}{Create diagnostic plots}
#'   \item{print.pgprobit}{Print model summaries}
#' }
#'
#' @references
#' Ding et al. (2024) \emph{Scalable Bayesian p-generalized probit and logistic
#' regression}. Advances in Data Analysis and Classification.
#' <doi:10.1007/s11634-024-00599-1>
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/zeyudsai/BayesPprobit}
#'   \item Report bugs at \url{https://github.com/zeyudsai/BayesPprobit/issues}
#' }
#'
#' @author
#' \strong{Maintainer}: Zeyu Ding \email{zeyu.ding@tu-dortmund.de}
#'
#' Authors:
#' \itemize{
#'   \item Katja Ickstadt
#'   \item Simon Omlor
#'   \item Alexander Munteanu
#' }
#'
#' @import stats
#' @importFrom graphics par plot abline
#' @importFrom methods is
#' @importFrom utils head tail str setTxtProgressBar txtProgressBar
#' @importFrom coda mcmc mcmc.list gelman.diag
#' @importFrom mvtnorm rmvnorm
#' @importFrom gnorm dgnorm pgnorm qgnorm
#' @importFrom Matrix crossprod solve
#' @importFrom sirt lq_fit
#'
#' @docType package
#' @name BayesPprobit-package
#' @aliases BayesPprobit
NULL
