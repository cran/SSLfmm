#' Initialize Parameters for a FMM from Labeled Subset
#'
#' Builds initial estimates \eqn{(\pi, \mu, \Sigma)} for a g-component Gaussian
#' mixture using only rows with observed labels in \code{zm}. Supports either a
#' shared covariance (\code{ncov = 1}) or class-specific covariances
#' (\code{ncov = 2}).
#'
#' @param dat A numeric matrix or data frame of features (n x p).
#' @param zm Integer vector of length n with class labels in \code{1:g}; use
#'   \code{NA} for unlabeled rows. Only labeled rows contribute to the
#'   initialization.
#' @param g Integer, number of mixture components.
#' @param ncov Integer, \code{1} for a shared covariance matrix, \code{2} for
#'   class-specific covariance matrices. Default \code{2}.
#' @param ridge Numeric, small diagonal ridge added to covariance(s) for
#'   numerical stability. Default \code{1e-6}.
#'
#' @return A list with
#' \itemize{
#'   \item \code{pi}: length-\code{g} vector of mixing proportions (summing to 1).
#'   \item \code{mu}: \code{p x g} matrix of class means (column \code{i} is \eqn{\mu_i}).
#'   \item \code{sigma}: if \code{ncov = 1}, a \code{p x p} shared covariance matrix;
#'         if \code{ncov = 2}, a \code{p x p x g} array of class-specific covariances.
#' }
#'
#' @details
#' If a class has zero or one labeled sample, its covariance is set to the global
#' empirical covariance (from labeled data) with a small ridge. Class means for
#' empty classes default to the global mean with a small jitter.
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 3; g <- 2
#' X <- matrix(rnorm(n*p), n, p)
#' z <- sample(c(1:g, NA), n, replace = TRUE, prob = c(0.4, 0.4, 0.2))
#' init <- initialestimate(X, z, g, ncov = 2)
#' str(init)
#'
#' @export
initialestimate <- function(dat, zm, g, ncov = 2, ridge = 1e-6) {
  # --- inputs ---
  if (!is.matrix(dat)) dat <- as.matrix(dat)
  if (!is.numeric(dat)) stop("`dat` must be numeric.")
  if (length(zm) != nrow(dat)) stop("`zm` must have length nrow(dat).")
  if (!ncov %in% c(1L, 2L)) stop("`ncov` must be 1 (shared) or 2 (class-specific).")
  if (g < 1L) stop("`g` must be >= 1.")
  
  # --- labeled subset ---
  labeled <- !is.na(zm)
  if (!any(labeled)) stop("No labeled rows in `zm`; cannot initialize.")
  Y <- dat[labeled, , drop = FALSE]
  k <- zm[labeled]
  p <- ncol(Y)
  n_lab <- length(k)
  
  # --- containers ---
  nn <- numeric(g)
  pi <- numeric(g)
  mu <- matrix(0, nrow = p, ncol = g)
  sigmaa <- array(0, dim = c(p, p, g))
  
  # --- global stats (fallbacks) ---
  glob_mean <- colMeans(Y)
  # use unbiased covariance (n-1) and ensure PSD with ridge
  glob_cov <- stats::cov(Y)
  if (any(!is.finite(glob_cov))) glob_cov <- diag(p)
  glob_cov <- glob_cov + diag(ridge, p)
  
  # --- per-class stats ---
  for (i in seq_len(g)) {
    idx <- which(k == i)
    nn[i] <- length(idx)
    pi[i] <- nn[i] / n_lab
    
    if (nn[i] >= 1) {
      Yi <- Y[idx, , drop = FALSE]
      mu[, i] <- colMeans(Yi)
      if (nn[i] >= 2) {
        # class covariance (unbiased)
        Si <- stats::cov(Yi)
        if (any(!is.finite(Si))) Si <- glob_cov
        sigmaa[, , i] <- Si + diag(ridge, p)
      } else {
        # only one sample -> fallback to global covariance
        sigmaa[, , i] <- glob_cov
      }
    } else {
      # empty class: mean fallback + small jitter; covariance fallback
      mu[, i] <- glob_mean + rnorm(p, sd = 1e-3)
      sigmaa[, , i] <- glob_cov
    }
  }
  
  # --- shared vs. class-specific covariance ---
  if (ncov == 1L) {
    # pooled within-class covariance (weighted by class counts)
    # sum over classes of (nn_i - 1) * S_i, then divide by (n_lab - g)
    # fall back to simple average if denominator < p
    denom <- sum(pmax(nn - 1, 0))
    if (denom > 0) {
      pooled <- Reduce(`+`, lapply(seq_len(g), function(i) pmax(nn[i] - 1, 0) * sigmaa[, , i])) / denom
    } else {
      pooled <- apply(sigmaa, c(1, 2), mean)
    }
    sigma <- pooled + diag(ridge, p)
  } else {
    # class-specific
    sigma <- sigmaa
  }
  
  # --- normalize pi (just in case) ---
  sumpi <- sum(pi)
  if (sumpi <= 0) {
    pi <- rep(1 / g, g)
  } else {
    pi <- pi / sumpi
  }
  
  list(pi = pi, mu = mu, sigma = sigma)
}
