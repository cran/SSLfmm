#' Posterior Cluster Probabilities for a Gaussian Mixture
#'
#' @description
#' Compute posterior membership probabilities \eqn{\tau_{ij}} for each
#' observation under a Gaussian mixture with either a shared covariance
#' matrix or component-specific covariances.
#'
#' @param dat Numeric matrix \eqn{n \times p}. Data matrix (rows = observations).
#' @param n Integer. Number of observations (checked against \code{nrow(dat)}).
#' @param p Integer. Number of variables (checked against \code{ncol(dat)}).
#' @param g Integer. Number of components.
#' @param pi Optional numeric vector length \eqn{g}. Mixing proportions (sum to 1).
#' @param mu Optional numeric matrix \eqn{p \times g}. Column \eqn{i} is mean of component \eqn{i}.
#' @param sigma Optional numeric matrix \eqn{p \times p} (shared covariance) or
#'   array \eqn{p \times p \times g} (component-specific covariances).
#' @param paralist Optional list with elements \code{pi}, \code{mu}, \code{sigma}.
#'   If provided, these take precedence over the corresponding explicit args.
#'
#' @return Numeric matrix \eqn{n \times g} of posterior probabilities \eqn{\tau_{ij}}.
#'
#' @examples
#'   n <- 100; p <- 2; g <- 2
#'   X <- matrix(rnorm(n*p), n, p)
#'   pi <- c(0.6, 0.4)
#'   mu <- cbind(c(0,0), c(1,1))
#'   Sig <- array(0, dim = c(p,p,g)); Sig[,,1] <- diag(p); Sig[,,2] <- diag(p)
#'   tau <- get_clusterprobs(X, n, p, g, pi = pi, mu = mu, sigma = Sig)
#'   head(tau)
#'
#' @importFrom mvtnorm dmvnorm
#' @export
get_clusterprobs <- function(dat, n, p, g, pi = NULL, mu = NULL, sigma = NULL, paralist = NULL){
  # Pull parameters from paralist if present
  if (!is.null(paralist)) {
    if (!is.null(paralist$pi))    pi    <- paralist$pi
    if (!is.null(paralist$mu))    mu    <- paralist$mu
    if (!is.null(paralist$sigma)) sigma <- paralist$sigma
  } else {
    paralist <- list(pi = pi, mu = mu, sigma = sigma)
  }

  # --- light checks ---
  if (!is.matrix(dat)) stop("dat must be a numeric matrix of size n x p.")
  if (nrow(dat) != n) stop("n does not match nrow(dat).")
  if (ncol(dat) != p) stop("p does not match ncol(dat).")
  if (is.null(pi) || is.null(mu) || is.null(sigma)) {
    stop("pi, mu, and sigma must be supplied (directly or via paralist).")
  }
  if (length(pi) != g) stop("length(pi) must equal g.")
  if (ncol(mu) != g || nrow(mu) != p) stop("mu must be a p x g matrix.")
  if (any(pi < 0) || abs(sum(pi) - 1) > 1e-8) stop("pi must be nonnegative and sum to 1.")

  # Determine covariance structure: shared (p x p) or per-component (p x p x g)
  sigma_dim <- dim(sigma)
  shared_cov <- FALSE
  if (length(sigma_dim) == 2) {
    if (!all(sigma_dim == c(p, p))) stop("sigma must be p x p if shared.")
    shared_cov <- TRUE
  } else if (length(sigma_dim) == 3) {
    if (!all(sigma_dim == c(p, p, g))) stop("sigma must be p x p x g if component-specific.")
  } else {
    stop("sigma must be either p x p (shared) or p x p x g (component-specific).")
  }

  # --- log densities for each component ---
  logdens <- matrix(NA_real_, nrow = n, ncol = g)
  if (shared_cov) {
    Sig <- as.matrix(sigma)
    for (k in seq_len(g)) {
      logdens[, k] <- mvtnorm::dmvnorm(dat, mean = mu[, k], sigma = Sig, log = TRUE)
    }
  } else {
    for (k in seq_len(g)) {
      logdens[, k] <- mvtnorm::dmvnorm(dat, mean = mu[, k], sigma = as.matrix(sigma[,,k]), log = TRUE)
    }
  }

  # Add log priors
  logprobs <- sweep(logdens, 2, log(pi), FUN = "+")

  # --- row-wise log-sum-exp and softmax normalization ---
  max_row <- apply(logprobs, 1, max)
  # subtract max per row (broadcast)
  centered <- logprobs - max_row
  lse <- max_row + log(rowSums(exp(centered)))
  # posterior probs
  tau <- exp(logprobs - lse)

  # Guard against tiny numerical drift: re-normalize rows to sum to 1
  tau <- tau / rowSums(tau)

  return(tau)
}
