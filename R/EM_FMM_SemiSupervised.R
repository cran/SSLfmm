#' EM for Semi-Supervised FMM with a Mixed-Missingness Mechanism (MCAR + entropy-based MAR)
#'
#' Runs an EM-like procedure that models a mixed-missingness mechanism:
#' unlabeled indicator \eqn{m_j} follows a mixed-missingness mechanism with MCAR (prob \eqn{\alpha})
#' and entropy-based MAR via a logistic link \eqn{q_j = \text{logit}^{-1}(\xi_0 + \xi_1 \log e_j)}.
#' Supports shared (\code{ncov = 1}) or class-specific (\code{ncov = 2}) covariance.
#'
#' @param data A data.frame or matrix with \eqn{p+2} columns:
#'   first \eqn{p} are features, then \code{missing} (0 = labelled, 1 = unlabelled),
#'   and \code{z} (class label for labelled rows; ignored otherwise).
#' @param g Integer, number of mixture components (classes).
#' @param init_res A list with initial parameters:
#'   \itemize{
#'     \item \code{pi}: numeric length-\code{g} (mixture weights, sum to 1)
#'     \item \code{mu}: list of length \code{g}, each length-\code{p} mean vector
#'     \item \code{Sigma}: if \code{ncov = 1}, a \code{p x p} matrix;
#'           if \code{ncov = 2}, a list of \code{g} \code{p x p} matrices
#'     \item \code{alpha}: scalar in (0,1)
#'     \item \code{xi}: numeric length-2, logistic coefficients \code{(xi0, xi1)}
#'   }
#' @param max_iter Integer, max EM iterations.
#' @param tol Convergence tolerance on log-likelihood increase.
#' @param ncov Integer covariance structure: \code{1} = shared/equal,
#'   \code{2} = class-specific/unequal.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed
#'   using \code{message()}. Default is \code{FALSE}.
#'
#' @details
#' This function expects the following helpers to be available:
#' \itemize{
#'   \item \code{pack_theta(pi_k, mu_k, Sigma_k, g, p, ncov)}
#'   \item \code{unpack_theta(theta, g, p, ncov)}
#'   \item \code{neg_loglik(theta, Y_all, m_j, Z_all, d2_yj, xi, alpha_k, unpacker)}
#'   \item \code{get_entropy(dat, n, p, g, paralist)} returning per-observation
#'         entropy-like values
#' }
#'
#' @return A list with elements:
#' \code{pi}, \code{mu}, \code{Sigma}, \code{xi}, \code{alpha}, \code{loglik},
#' and \code{ncov}.
#'
#' @examples
#' \donttest{
#'   ## Toy example using a simple partially-labelled dataset
#'   set.seed(1)
#'
#'   ## 1) Construct an n x (p+2) partially-labelled dataset:
#'   ##    first p columns = features, then 'missing' and 'z'
#'   n <- 100; p <- 2; g <- 2
#'
#'   X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#'
#'   ## missing: 0 = labelled, 1 = unlabelled
#'   missing <- rbinom(n, size = 1, prob = 0.3)
#'
#'   ## z: observed class labels for labelled rows, NA for unlabelled
#'   z <- rep(NA_integer_, n)
#'   z[missing == 0] <- sample(1:g, sum(missing == 0), replace = TRUE)
#'
#'   sim_dat <- data.frame(X, missing = missing, z = z)
#'
#'   ## 2) Warm-up initialisation using the complete-data initializer
#'   init <- EM_FMM_SemiSupervised_Complete_Initial(
#'     data = sim_dat,
#'     g    = g,
#'     ncov = 1
#'   )
#'
#'   ## 3) Run the main EM algorithm (small number of iterations)
#'   fit <- EM_FMM_SemiSupervised(
#'     data     = sim_dat,
#'     g        = g,
#'     init_res = init,
#'     ncov     = 1,
#'     max_iter = 5,
#'     verbose  = FALSE
#'   )
#'
#'   str(fit)
#' }
#'
#' @export
EM_FMM_SemiSupervised <- function(data,
                                  g        = 2,
                                  init_res,
                                  max_iter = 5,
                                  tol      = 1e-6,
                                  ncov     = 1,
                                  verbose  = FALSE) {
  if (!ncov %in% c(1, 2)) stop("ncov must be 1 (shared) or 2 (class-specific).")
  
  # pull in needed namespaces if not already attached
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("Package 'mvtnorm' is required.")
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("Package 'matrixStats' is required.")
  
  data <- as.data.frame(data)
  
  # assume last two columns are 'missing' and 'z'
  p <- ncol(data) - 2L
  Y_labelled   <- data[data[["missing"]] == 0, seq_len(p), drop = FALSE]
  Y_unlabelled <- data[data[["missing"]] != 0, seq_len(p), drop = FALSE]
  Z_labelled   <- data[data[["missing"]] == 0, "z"]
  
  Y_all     <- rbind(Y_labelled, Y_unlabelled)
  n_label   <- nrow(Y_labelled)
  n_unlabel <- nrow(Y_unlabelled)
  n         <- n_label + n_unlabel
  p         <- ncol(Y_labelled)
  
  m_j   <- c(rep(0L, n_label), rep(1L, n_unlabel))
  Z_all <- c(Z_labelled, rep(NA_integer_, n_unlabel))
  
  eps        <- 1e-10
  loglik_old <- -Inf
  loglik1    <- NA_real_
  
  # ---- initial params ----
  pi_k    <- init_res$pi
  mu_k    <- init_res$mu
  Sigma_k <- init_res$Sigma
  alpha_k <- init_res$alpha
  xi      <- init_res$xi
  
  # helper: build p x p x g array from Sigma (replicate if shared)
  make_sigma_array <- function(Sigma_k, g) {
    if (is.list(Sigma_k)) {
      simplify2array(Sigma_k)
    } else {
      A <- array(0, dim = c(p, p, g))
      for (k in seq_len(g)) A[, , k] <- Sigma_k
      A
    }
  }
  
  for (iter in seq_len(max_iter)) {
    # ----- E-step -----
    Sigma_arr <- make_sigma_array(Sigma_k, g)
    ests <- list(
      pi    = pi_k,
      mu    = do.call(cbind, mu_k),  # p x g
      sigma = Sigma_arr              # p x p x g
    )
    
    d2_yj <- log(get_entropy(dat = Y_all, n = n, p = p, g = g, paralist = ests) + eps)
    
    q_j   <- plogis(xi[1L] + xi[2L] * d2_yj)
    denom <- alpha_k + (1 - alpha_k) * q_j
    m1j   <- alpha_k / denom
    m2j   <- (1 - alpha_k) * q_j / denom
    m1j_k <- ifelse(m_j == 1L, m1j, 0)
    m2j_k <- ifelse(m_j == 1L, m2j, 0)
    
    # ----- M-step: alpha, xi -----
    alpha_k <- mean(m1j_k)
    glm_df  <- data.frame(m2j_k = m2j_k, d2_yj = d2_yj)
    xi <- withCallingHandlers(
      {
        stats::coef(stats::glm(
          m2j_k ~ d2_yj,
          data    = glm_df,
          family  = stats::quasibinomial(link = "logit"),
          control = stats::glm.control(maxit = 100)
        ))
      },
      warning = function(w) invokeRestart("muffleWarning")
    )
    
    # ----- M-step: pi, mu, Sigma via optimization -----
    theta0 <- pack_theta(pi_k, mu_k, Sigma_k, g, p, ncov = ncov)
    
    fit <- stats::nlminb(
      start     = theta0,
      objective = function(th) {
        neg_loglik(th, Y_all, m_j, Z_all,
                   d2_yj, xi, alpha_k,
                   function(x) unpack_theta(x, g, p, ncov = ncov))
      }
    )
    
    loglik1 <- -fit$objective
    newp    <- unpack_theta(fit$par, g, p, ncov = ncov)
    
    pi_k    <- newp$pi
    mu_k    <- newp$mu
    Sigma_k <- newp$sigma
    
    if (verbose) {
      msg <- sprintf(
        "Iter %d: full_nll = %.6f | alpha = %.4f | xi0 = %.4f | xi1 = %.4f | ncov = %d",
        iter, loglik1, alpha_k, xi[1L], xi[2L], ncov
      )
      message(msg)
    }
    
    if (abs(loglik1 - loglik_old) < tol) {
      if (verbose) message("Converged.")
      break
    }
    loglik_old <- loglik1
  }
  
  list(
    pi     = pi_k,
    mu     = mu_k,
    Sigma  = Sigma_k,
    xi     = xi,
    alpha  = alpha_k,
    loglik = loglik1,
    ncov   = ncov
  )
}
