#' @title Package Imports for BayesGProbit
#' @description Centralized management of package imports
#' @name imports
#' @keywords internal
NULL

#----------------------
# Statistics Packages
#----------------------

#' @importFrom stats
#'   runif
#'   rnorm
#'   rbinom
#'   rgamma
#'   pgamma
#'   qgamma
#'   density
#'   quantile
#'   sd
#'   var
#'   cov
#'   cor
#'   median
#'   IQR
#'   pnorm
#'   qnorm
#'   integrate
#'   optimize
NULL

#----------------------
# Matrix Computations
#----------------------

#' @importFrom mvtnorm
#'   rmvnorm
#'   dmvnorm
NULL

#' @importFrom Matrix
#'   crossprod
#'   tcrossprod
#'   qr
#'   qr.Q
#'   qr.R
#'   solve
NULL

#----------------------
# Graphics & Plotting
#----------------------

#' @importFrom graphics
#'   par
#'   plot
#'   points
#'   lines
#'   abline
#'   legend
#'   title
#'   axis
#'   box
#'   hist
#'   mtext
NULL

#' @importFrom grDevices
#'   dev.off
#'   pdf
#'   png
#'   rgb
#'   col2rgb
#'   colorRampPalette
NULL

#----------------------
# Methods & Classes
#----------------------

#' @importFrom methods
#'   is
#'   new
#'   slot
#'   slotNames
#'   validObject
NULL

#----------------------
# MCMC Diagnostics
#----------------------

#' @importFrom coda
#'   mcmc
#'   mcmc.list
#'   gelman.diag
#'   effectiveSize
#'   geweke.diag
NULL

#' @importFrom truncnorm
#'   rtruncnorm
#'   dtruncnorm
#'   ptruncnorm
#'   qtruncnorm
NULL

#----------------------
# Data Manipulation
#----------------------

#' @importFrom utils
#'   head
#'   tail
#'   str
#'   write.csv
#'   read.csv
#'   setTxtProgressBar
#'   txtProgressBar
NULL

#----------------------
# Optional Packages
#----------------------

# These are in Suggests, not Imports
if(getRversion() >= "3.6.0") {
  #' @importFrom ggplot2
  #'   ggplot
  #'   aes
  #'   geom_point
  #'   geom_line
  #'   theme_bw
  NULL
}

#----------------------
# Package Documentation
#----------------------

#' @import stats
"_PACKAGE"

# Make internal functions available within package
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Any package initialization code
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("BayesGProbit version ", utils::packageVersion("BayesGProbit"))
}

# Register S3 methods
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Register S3 methods
  s3_register("graphics::plot", "pgprobit")
  s3_register("base::print", "pgprobit")
  s3_register("base::summary", "pgprobit")
}

# Utility function to safely register S3 methods
s3_register <- function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  pieces <- strsplit(generic, "::")[[1]]
  if(length(pieces) == 1) {
    package <- topenv(parent.frame())$`__package__`
    pieces <- c(package, pieces)
  }

  if(is.null(method)) {
    method <- get(paste0(pieces[2], ".", class), envir = parent.frame())
  }

  registerS3method(pieces[2], class, method, pieces[1])
}
