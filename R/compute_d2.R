#' Squared Discriminant Score for Two-Group LDA (Equal Covariance)
#'
#' Computes \eqn{d^2} where \eqn{d = \beta_0 + y^\top \beta_1} for two classes
#' under the common-covariance (LDA) model:
#' \deqn{\beta_1 = \Sigma^{-1}(\mu_1 - \mu_2), \quad
#'       \beta_0 = \log(\pi_1/\pi_2) - \tfrac{1}{2}(\mu_1 + \mu_2)^\top \Sigma^{-1}(\mu_1 - \mu_2).}
#'
#' @param y Numeric vector (length \eqn{p}) or numeric matrix with \eqn{p} columns;
#'   rows are observations.
#' @param mu1,mu2 Numeric vectors of length \eqn{p}: class means.
#' @param Sigma_inv \eqn{p \times p} numeric precision matrix (inverse covariance).
#' @param pi1,pi2 Positive scalars: class prior probabilities (need not sum to 1).
#'
#' @return If \code{y} is a vector, a single numeric \eqn{d^2}. If \code{y} is a matrix,
#'   a numeric vector of \eqn{d^2} values for each row.
#'
#' @examples
#' set.seed(1)
#' mu1 <- c(0, 0); mu2 <- c(1, 1)
#' S <- matrix(c(1, .2, .2, 1), 2, 2)
#' Sigma_inv <- solve(S)
#' x <- c(0.5, -0.2)
#' compute_d2(x, mu1, mu2, Sigma_inv, pi1 = 0.6, pi2 = 0.4)
#'
#' X <- matrix(rnorm(10 * 2), ncol = 2)
#' compute_d2(X, mu1, mu2, Sigma_inv, 0.5, 0.5)
#'
#' @export
compute_d2 <- function(y, mu1, mu2, Sigma_inv, pi1, pi2) {
  # Basic checks
  if (length(mu1) != length(mu2)) stop("mu1 and mu2 must have same length.")
  p <- length(mu1)
  if (!all(dim(Sigma_inv) == c(p, p))) stop("Sigma_inv must be p x p.")
  if (any(c(pi1, pi2) <= 0)) stop("pi1 and pi2 must be positive.")
  
  # Ensure y is a matrix with p columns
  if (is.vector(y)) y <- matrix(y, nrow = 1)
  if (ncol(y) != p) stop("Number of columns of y must equal length(mu1).")
  
  # Discriminant parameters
  beta1 <- Sigma_inv %*% (mu1 - mu2)
  beta0 <- as.numeric(log(pi1 / pi2) - 0.5 * t(mu1 + mu2) %*% beta1)
  
  # d for each row of y, then square
  d <- as.numeric(y %*% beta1) + beta0
  d^2
}
