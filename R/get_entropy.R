#' Per-Row Entropy of Posterior Cluster Probabilities
#'
#' @description
#' Compute the Shannon entropy (in nats) of posterior membership probabilities
#' for each observation, given a Gaussian mixture model. Posterior probabilities
#' \eqn{\tau_{ij}} are obtained via \code{get_clusterprobs()}.
#'
#' @param dat Numeric matrix \eqn{n \times p}. The data matrix.
#' @param n Integer. Number of rows in \code{dat}.
#' @param p Integer. Number of columns (features) in \code{dat}.
#' @param g Integer. Number of mixture components.
#' @param pi Optional numeric vector length \eqn{g}. Mixing proportions (sum to 1).
#' @param mu Optional numeric matrix \eqn{p \times g}. Column \eqn{j} is the mean of component \eqn{j}.
#' @param sigma Optional numeric matrix \eqn{p \times p} (shared covariance) or array
#'   \eqn{p \times p \times g} (component-specific covariances).
#' @param paralist Optional list with elements \code{pi}, \code{mu}, \code{sigma}.
#'   If supplied, these take precedence over the corresponding explicit arguments.
#'
#' @return Numeric vector of length \eqn{n}, the entropy per observation (in nats).
#'
#' @examples
#'   # Suppose you have get_clusterprobs(), and a fitted/pi, mu, sigma:
#'   n <- 100; p <- 2; g <- 2
#'   X <- matrix(rnorm(n * p), n, p)
#'   pi <- c(0.6, 0.4)
#'   mu <- cbind(c(0,0), c(1,1))
#'   Sigma <- array(0, dim = c(p, p, g))
#'   Sigma[,,1] <- diag(p); Sigma[,,2] <- diag(p)
#'   ent <- get_entropy(dat = X, n = n, p = p, g = g, pi = pi, mu = mu, sigma = Sigma)
#'   summary(ent)
#'
#' @seealso \code{\link{get_clusterprobs}}
#' @export
get_entropy <- function(dat, n, p, g, pi = NULL, mu = NULL, sigma = NULL, paralist = NULL) {
  # --- pull parameters from paralist if provided ---
  if (!is.null(paralist)) {
    if (!is.null(paralist$pi))    pi    <- paralist$pi
    if (!is.null(paralist$mu))    mu    <- paralist$mu
    if (!is.null(paralist$sigma)) sigma <- paralist$sigma
  }
  dat<-as.matrix(dat)
  # --- light checks (avoid heavy validation to keep this fast) ---
  if (!is.matrix(dat)) stop("dat must be a numeric matrix of size n x p.")
  if (nrow(dat) != n) stop("n does not match nrow(dat).")
  if (ncol(dat) != p) stop("p does not match ncol(dat).")
  if (is.null(pi) || is.null(mu) || is.null(sigma)) {
    stop("pi, mu, and sigma must be supplied (either directly or via paralist).")
  }
  
  # --- posterior probabilities (n x g) ---
  tau <- get_clusterprobs(dat = dat, n = n, p = p, g = g, mu = mu, sigma = sigma, pi = pi)
  
  # numerical stability: clamp tau to at least eps, preserves 0*log(0)=0 behavior
  eps <- .Machine$double.eps
  tau_safe <- pmax(tau, eps)
  entropy <- -rowSums(tau * log(tau_safe))
  
  return(entropy)
}
