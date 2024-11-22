#' @title Construct Coreset for p-Generalized Probit Models
#' @description Constructs coresets for efficient data compression using various sampling methods
#' @param X Design matrix
#' @param y Binary response vector
#' @param coreset_size Desired size of coreset
#' @param method Sampling method: "leverage" (default), "uniform", or "oneshot"
#' @param sketch_size Size of sketching matrix (default: ncol(X)^2)
#' @param p Shape parameter for p-norm (default: 2)
#' @return A list containing:
#'   \item{indices}{Selected data indices}
#'   \item{weights}{Importance weights}
#'   \item{scores}{Sampling scores}
#' @export
compute_coreset <- function(X, y, coreset_size,
                            method = c("leverage", "uniform", "oneshot"),
                            sketch_size = NULL,
                            p = 2) {

  # Input validation
  if(!methods::is(X, "matrix")) stop("X must be a matrix")
  if(!is.numeric(y) || !all(y %in% c(0,1)))
    stop("y must be binary (0/1)")
  if(length(y) != nrow(X))
    stop("Length of y must match rows of X")

  method <- match.arg(method)

  n <- nrow(X)
  d <- ncol(X)

  # Default sketch size
  if(is.null(sketch_size)) {
    sketch_size <- d^2
  }

  # Compute sampling scores based on method
  scores <- switch(method,
                   "leverage" = compute_leverage_scores(X, p, sketch_size),
                   "uniform" = rep(1/n, n),
                   "oneshot" = compute_oneshot_scores(X, p, sketch_size)
  )

  # Sample indices
  indices <- sample.int(n,
                        size = coreset_size,
                        prob = scores,
                        replace = FALSE)

  # Compute importance weights
  weights <- 1/(coreset_size * scores[indices])

  list(
    indices = indices,
    weights = weights,
    scores = scores
  )
}

#' @title Compute Leverage Scores
#' @description Computes p-norm leverage scores using sketching
#' @keywords internal
compute_leverage_scores <- function(X, p, sketch_size) {
  n <- nrow(X)
  d <- ncol(X)

  # Compute sketch
  X_sketch <- compute_sketch(X, sketch_size, p)

  # QR decomposition of sketch
  R <- Matrix::qr.R(Matrix::qr(X_sketch))

  # Compute scores
  XR <- X %*% Matrix::solve(R)
  scores <- numeric(n)

  for(i in 1:n) {
    scores[i] <- norm(XR[i,], type = "2")^p
  }

  # Add baseline probability
  scores <- scores + 1/n

  # Normalize
  scores/sum(scores)
}

#' @title Compute One-shot Coreset Scores
#' @description Computes sampling scores valid for multiple p values
#' @keywords internal
compute_oneshot_scores <- function(X, p_max, sketch_size) {
  n <- nrow(X)
  delta <- 1/log(n)

  # Grid of p values
  p_grid <- seq_geom(1, p_max, by = 1 + delta)

  # Compute scores for each p
  scores_all <- matrix(0, nrow = n, ncol = length(p_grid))
  for(i in seq_along(p_grid)) {
    scores_all[,i] <- compute_leverage_scores(X, p_grid[i], sketch_size)
  }

  # Take maximum over all p values
  scores <- apply(scores_all, 1, max)

  # Normalize
  scores/sum(scores)
}

#' @title Compute Matrix Sketch
#' @description Computes Count-Sketch or Gaussian sketch of input matrix
#' @keywords internal
compute_sketch <- function(X, sketch_size, p) {
  n <- nrow(X)
  d <- ncol(X)

  # Generate random indices and signs
  h <- sample.int(sketch_size, n, replace = TRUE)
  s <- sample(c(-1,1), n, replace = TRUE)

  # For p != 2, adjust signs
  if(p != 2) {
    lambda <- stats::rexp(n, 1)
    s <- s/(lambda^(1/p))
  }

  # Construct sketch matrix
  S <- Matrix::sparseMatrix(
    i = h,
    j = 1:n,
    x = s,
    dims = c(sketch_size, n)
  )

  # Return sketched matrix
  as.matrix(S %*% X)
}

#' @title Generate Geometric Sequence
#' @description Generates sequence with geometric spacing
#' @keywords internal
seq_geom <- function(from, to, by) {
  exp(seq(log(from), log(to), log(by)))
}

