#' Bayes' Rule Classifier
#'
#' Classifier specified by Bayes' rule.
#' Assigns \eqn{\arg\max_i \{ \log(\pi_i) + \log \mathcal{N} (y_j \mid \mu_i, \Sigma_i) \}}.
#'
#' @param dat An \eqn{n \times p} numeric matrix (or a length-\eqn{p} numeric vector,
#'   treated as \eqn{1 \times p}).
#' @param p Integer; dimension of the observation vector.
#' @param g Integer; number of Gaussian components/classes.
#' @param pi Numeric length-\eqn{g} vector of mixing proportions (must sum to 1).
#' @param mu Either:
#'   \itemize{
#'     \item A \eqn{p \times g} numeric matrix (column \eqn{i} is \eqn{\mu_i}), or
#'     \item A length-\eqn{g} list of length-\eqn{p} numeric vectors.
#'   }
#' @param sigma Either:
#'   \itemize{
#'     \item A \eqn{p \times p} covariance matrix (shared), or
#'     \item A length-\eqn{g} list of \eqn{p \times p} covariance matrices, or
#'     \item A \eqn{p \times p \times g} numeric array.
#'   }
#' @param paralist Optional list with elements \code{pi}, \code{mu}, \code{sigma}
#'   (overrides explicit args if provided and non-\code{NULL}).
#'
#' @return An integer vector of length \eqn{n} with predicted class labels in \eqn{1:g}.
#'
#' @examples
#' # Minimal example with list-style mu and sigma:
#' set.seed(1)
#' p <- 2; g <- 2
#' pi <- c(0.5, 0.5)
#' mu <- list(`1`=c(0.96, 0.02), `2`=c(-1.02, -0.03))
#' sigma <- list(
#'   matrix(c(0.9417379,0.5447264, 0.5447264,0.9811853), 2, 2, byrow=TRUE),
#'   matrix(c(0.9984812,0.3314474, 0.3314474,1.1316865), 2, 2, byrow=TRUE)
#' )
#' X <- mvtnorm::rmvnorm(5, mean = mu[[1]], sigma = sigma[[1]])
#' bayesclassifier(X, p=p, g=g, pi=pi, mu=mu, sigma=sigma)
#'
#' @importFrom mvtnorm dmvnorm
#' @export
bayesclassifier <- function(dat, p, g, pi = NULL, mu = NULL, sigma = NULL, paralist = NULL) {
  
  # ---- pull from paralist if provided (only override non-NULL fields) ----
  if (!is.null(paralist)) {
    if (!is.null(paralist$pi))    pi    <- paralist$pi
    if (!is.null(paralist$mu))    mu    <- paralist$mu
    if (!is.null(paralist$sigma)) sigma <- paralist$sigma
  } else {
    paralist <- list(pi = pi, mu = mu, sigma = sigma)
  }
  
  # ---- coerce dat to \eqn{n \times p} matrix; allow vector input ----
  if (is.null(dim(dat))) {
    dat <- matrix(dat, nrow = 1)
  } else {
    dat <- as.matrix(dat)
  }
  n <- nrow(dat)
  if (ncol(dat) != p) stop("ncol(dat) must equal p.")
  
  # ---- basic checks on pi ----
  if (is.null(pi) || length(pi) != g)
    stop("pi must be provided and have length g.")
  if (any(pi < 0) || abs(sum(pi) - 1) > 1e-8)
    stop("pi must be nonnegative and sum to 1.")
  
  # ---- coerce mu: list -> p times g matrix ----
  if (is.null(mu)) stop("mu must be provided.")
  if (is.list(mu)) {
    if (length(mu) != g) stop("mu list length must be g.")
    mu <- do.call(cbind, mu)
  }
  mu <- as.matrix(mu)
  if (!all(dim(mu) == c(p, g)))
    stop("mu must be a p times g matrix (after coercion).")
  
  # ---- coerce sigma: list -> array p times p times g; or allow shared p times p ----
  if (is.null(sigma)) stop("sigma must be provided.")
  if (is.list(sigma)) {
    # list of p×p -> array
    if (length(sigma) != g) stop("sigma list length must be g.")
    pchk <- nrow(sigma[[1]])
    if (pchk != p) stop("sigma list elements must be p times p.")
    arr <- array(NA_real_, dim = c(p, p, g))
    for (k in seq_len(g)) {
      if (!all(dim(sigma[[k]]) == c(p, p)))
        stop("each sigma[[k]] must be p times p.")
      arr[,,k] <- sigma[[k]]
    }
    sigma <- arr
  }
  
  # determine covariance structure
  sd <- dim(sigma)
  shared_cov <- FALSE
  if (length(sd) == 2) {
    if (!all(sd == c(p, p))) stop("sigma must be p times p if shared.")
    shared_cov <- TRUE
  } else if (length(sd) == 3) {
    if (!all(sd == c(p, p, g))) stop("sigma must be p times p times g if component-specific.")
  } else {
    stop("sigma must be either pxp (shared) or pxpxg (component-specific).")
  }
  
  # ---- compute log densities per component ----
  logdens <- matrix(NA_real_, nrow = n, ncol = g)
  if (shared_cov) {
    Sig <- as.matrix(sigma)
    for (k in seq_len(g)) {
      logdens[, k] <- mvtnorm::dmvnorm(dat, mean = mu[, k], sigma = Sig, log = TRUE)
    }
  } else {
    for (k in seq_len(g)) {
      logdens[, k] <- mvtnorm::dmvnorm(dat, mean = mu[, k], sigma = sigma[,,k], log = TRUE)
    }
  }
  
  # ---- add log priors and take argmax (no need to normalize) ----
  logprobs <- sweep(logdens, 2, log(pi), `+`)
  clust <- max.col(logprobs, ties.method = "first")
  
  return(clust)
}
