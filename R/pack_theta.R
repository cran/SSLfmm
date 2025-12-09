#' Pack FMM Parameters into a Vector
#'
#' Packs mixture weights, means, and covariance(s) into a single numeric vector.
#' Uses the last component as the baseline for mixture weights (g-1 logits stored).
#'
#' @param pi_k Numeric vector of length g with mixture weights (positive, sum to 1).
#' @param mu_k List of length g; each element a numeric vector of length p (component means).
#' @param Sigma Covariance: if \code{ncov = 1}, a single p x p matrix;
#'   if \code{ncov = 2}, a list of g p x p matrices.
#' @param g Integer: number of components.
#' @param p Integer: dimension.
#' @param ncov Integer: covariance structure; 1 for shared covariance, 2 for class-specific.
#'
#' @return Numeric vector with parameters packed.
#' @export
pack_theta <- function(pi_k, mu_k, Sigma, g, p, ncov = 1) {
  if (length(pi_k) != g) stop("pi_k must have length g.")
  if (!isTRUE(all.equal(sum(pi_k), 1, tolerance = 1e-8))) stop("sum(pi_k) must be 1.")
  if (length(mu_k) != g) stop("mu_k must be a list of length g.")
  if (!all(vapply(mu_k, length, 1L) == p)) stop("Each mu_k[[k]] must have length p.")
  if (!ncov %in% c(1, 2)) stop("ncov must be 1 (shared) or 2 (class-specific).")
  
  # logit transform for pi (last class baseline)
  logit_pi <- log(pi_k[-g] / pi_k[g])
  mu_vec <- unlist(mu_k, use.names = FALSE)
  
  vech <- function(M) M[lower.tri(M, diag = TRUE)]
  if (ncov == 1) {
    if (!is.matrix(Sigma) || any(dim(Sigma) != c(p, p))) stop("Sigma must be p x p matrix.")
    sigma_vec <- vech(Sigma)
    c(logit_pi, mu_vec, sigma_vec)
  } else {
    if (!is.list(Sigma) || length(Sigma) != g) stop("Sigma must be a list of length g.")
    sigma_vec <- unlist(lapply(Sigma, vech), use.names = FALSE)
    c(logit_pi, mu_vec, sigma_vec)
  }
}