#' Unpack FMM Parameter Vector
#'
#' Unpacks mixture weights, means, and covariance(s) from a parameter vector.
#'
#' @param theta Numeric vector as returned by \code{pack_theta()}.
#' @param g Integer: number of components.
#' @param p Integer: dimension.
#' @param ncov Integer: covariance structure; 1 for shared covariance, 2 for class-specific.
#'
#' @return A list with:
#' \describe{
#'   \item{pi}{Mixture weights (length g).}
#'   \item{mu}{List of g mean vectors.}
#'   \item{sigma}{Shared covariance matrix (ncov=1) or list of g covariance matrices (ncov=2).}
#' }
#' @export
unpack_theta <- function(theta, g, p, ncov = 1) {
  if (!ncov %in% c(1, 2)) stop("ncov must be 1 (shared) or 2 (class-specific).")
  
  idx <- 0L
  logit_pi <- theta[seq_len(g - 1L)]
  idx <- idx + (g - 1L)
  
  # Stable softmax with last class baseline
  max_l <- max(c(logit_pi, 0))
  exp_l <- exp(c(logit_pi, 0) - max_l)
  pi_k <- exp_l / sum(exp_l)
  
  mu_slice <- theta[(idx + 1L):(idx + g * p)]
  idx <- idx + g * p
  mu_k <- lapply(seq_len(g), function(k) {
    mu_slice[((k - 1L) * p + 1L):(k * p)]
  })
  names(mu_k) <- as.character(seq_len(g))
  
  ivech <- function(v, p) {
    M <- matrix(0, p, p)
    M[lower.tri(M, diag = TRUE)] <- v
    M + t(M) - diag(diag(M))
  }
  
  n_sigma <- p * (p + 1L) / 2L
  if (ncov == 1) {
    sigma_vals <- theta[(idx + 1L):(idx + n_sigma)]
    Sigma_shared <- ivech(sigma_vals, p)
    sigma_out <- Sigma_shared
  } else {
    sigma_vals <- theta[(idx + 1L):(idx + g * n_sigma)]
    sigma_out <- lapply(seq_len(g), function(k) {
      ivech(sigma_vals[((k - 1L) * n_sigma + 1L):(k * n_sigma)], p)
    })
  }
  
  list(pi = pi_k, mu = mu_k, sigma = sigma_out)
}