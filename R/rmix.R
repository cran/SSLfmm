#' Draw from a Gaussian Mixture Model
#'
#' @description
#' Generate i.i.d. samples from a finite Gaussian mixture with either a shared
#' covariance matrix or component-specific covariance matrices.
#'
#' @param n Integer. Number of observations to generate.
#' @param pi Numeric vector of length \eqn{g}. Mixing proportions (must sum to 1).
#' @param mu Numeric matrix \eqn{p \times g}. Column \eqn{j} is the mean for component \eqn{j}.
#' @param sigma Either a numeric matrix \eqn{p \times p} (shared covariance),
#'   or a numeric array \eqn{p \times p \times g} (component-specific covariances).
#' @param seed_number Integer. Seed for reproducibility.
#'
#' @return A list with:
#' \item{Y}{Numeric matrix \eqn{n \times p} of generated features.}
#' \item{Z}{Numeric matrix \eqn{n \times g} of one-hot component indicators.}
#' \item{clust}{Integer vector \eqn{n}, the component labels in \code{1:g}.}
#'
#' @examples
#'   set.seed(1)
#'   g  <- 2; p <- 2
#'   pi <- c(0.5, 0.5)
#'   mu <- cbind(c(1,0), c(-1,0))
#'   Sigma <- diag(p)
#'   out <- rmix(500, pi, mu, Sigma, seed_number = 123)
#'   str(out)
#'
#' @importFrom mvtnorm rmvnorm
#' @export
rmix <- function(n, pi, mu, sigma, seed_number) {
  # ---- checks ----
  if (!is.numeric(n) || length(n) != 1 || n < 0 || n != as.integer(n)) {
    stop("n must be a non-negative integer.")
  }
  if (!is.numeric(pi) || any(pi < 0) || abs(sum(pi) - 1) > 1e-8) {
    stop("pi must be nonnegative and sum to 1.")
  }
  if (!is.matrix(mu)) stop("mu must be a numeric matrix (p x g).")
  g <- length(pi)
  p <- nrow(mu)
  if (ncol(mu) != g) stop("ncol(mu) must equal length(pi).")
  
  # sigma: allow p x p (shared) or p x p x g (per-component)
  sigma_dim <- dim(sigma)
  if (length(sigma_dim) == 2) {
    if (!all(sigma_dim == c(p, p))) stop("sigma must be p x p if shared.")
    shared_cov <- TRUE
  } else if (length(sigma_dim) == 3) {
    if (!all(sigma_dim == c(p, p, g))) stop("sigma must be p x p x g.")
    shared_cov <- FALSE
  } else {
    stop("sigma must be either a p x p matrix or a p x p x g array.")
  }
  
  if (!is.numeric(seed_number) || length(seed_number) != 1) {
    stop("seed_number must be a single numeric value.")
  }
  
  # ---- sampling ----
  set.seed(seed_number)
  
  # component counts (length g; guaranteed to include all components)
  nn <- as.vector(rmultinom(1, n, prob = pi))
  
  # preallocate list for blocks
  X_blocks <- vector("list", g)
  for (j in seq_len(g)) {
    nj <- nn[j]
    if (nj > 0) {
      Sigma_j <- if (shared_cov) sigma else sigma[ , , j, drop = FALSE][,,1]
      X_blocks[[j]] <- mvtnorm::rmvnorm(nj, mean = mu[, j], sigma = Sigma_j)
    } else {
      X_blocks[[j]] <- matrix(numeric(0), nrow = 0, ncol = p)
    }
  }
  
  # stack and label
  X <- do.call(rbind, X_blocks)
  y <- rep(seq_len(g), times = nn)
  
  # shuffle
  if (n > 0) {
    rperm <- sample.int(n)
    X <- X[rperm, , drop = FALSE]
    y <- y[rperm]
  } else {
    X <- matrix(numeric(0), nrow = 0, ncol = p)
    y <- integer(0)
  }
  
  # one-hot (n x g)
  Z <- if (n > 0) {
    diag(g)[y, , drop = FALSE]
  } else {
    matrix(numeric(0), nrow = 0, ncol = g)
  }
  
  list(Y = X, Z = Z, clust = y)
}
