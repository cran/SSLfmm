#' Quick Initializer for alpha, xi, and Mixture Parameters
#'
#' @description
#' Provides rough initial estimates of the missingness parameters \code{alpha}
#' and \code{xi}, together with mixture parameters \code{pi}, \code{mu}, and
#' covariance matrices, using a lightweight EM-style routine. The covariance
#' structure is chosen automatically based on \code{Sigma_init}:
#'
#' \itemize{
#'   \item If \code{Sigma_init} is a \eqn{p \times p} matrix, a
#'         \strong{shared (equal)} covariance is used.
#'   \item If \code{Sigma_init} is a list of length \code{g} of
#'         \eqn{p \times p} matrices or a \eqn{p \times p \times g} array,
#'         \strong{class-specific (unequal)} covariances are used.
#'   \item If \code{Sigma_init} is \code{NULL}, a shared covariance is
#'         estimated from the labeled data.
#' }
#'
#' This function is intended as a fast, heuristic initializer rather than a
#' final estimator for the mixed missingness model.
#'
#' @param Y_labelled Numeric matrix of labeled observations (\eqn{n_L \times p}).
#' @param Z_labelled Integer vector of class labels in \code{1:g} for
#'   \code{Y_labelled}.
#' @param Y_unlabelled Numeric matrix of unlabeled observations
#'   (\eqn{n_U \times p}).
#' @param g Integer, number of mixture components (default \code{2}).
#' @param pi_init Optional numeric length-\code{g} vector of initial mixing
#'   proportions.
#' @param mu_init Optional list of length \code{g} of initial mean vectors
#'   (each of length \code{p}).
#' @param Sigma_init Optional initial covariance:
#'   a \eqn{p \times p} matrix (shared), or a list of \code{g}
#'   \eqn{p \times p} matrices, or a \eqn{p \times p \times g} array
#'   (class-specific).
#' @param alpha_init Numeric in \eqn{(0,1)}, initial MCAR proportion
#'   (default \code{0.01}).
#' @param warm_up_iter Integer, number of warm-up EM iterations used to
#'   refine the quick initial estimates (default \code{50}).
#' @param tol Convergence tolerance on \code{alpha} (default \code{1e-6}).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{pi} - length-\code{g} vector of mixing proportions.
#'   \item \code{mu} - list of \code{g} mean vectors.
#'   \item \code{Sigma} - shared \eqn{p \times p} matrix (equal-Sigma) or list
#'         of \code{g} matrices (unequal-Sigma).
#'   \item \code{xi} - length-2 numeric vector \code{c(xi0, xi1)} from the
#'         logistic MAR model.
#'   \item \code{alpha} - estimated MCAR proportion.
#'   \item \code{gamma} - \eqn{n \times g} responsibility matrix.
#'   \item \code{d2_yj} - numeric vector of entropy-based scores used in the
#'         missingness model.
#' }
#'
#' @note
#' This is a heuristic warm-up routine. It requires helper functions
#' \code{normalise_logprob()} and \code{get_entropy()} to be available in the
#' package namespace.
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats glm plogis coef cov
#' @export
EM_FMM_SemiSupervised_Initial <- function(Y_labelled, Z_labelled, Y_unlabelled, g = 2, 
                                          pi_init = NULL, mu_init = NULL, Sigma_init = NULL,
                                          alpha_init = 0.01,
                                          warm_up_iter = 50, tol = 1e-6) {
  
  n_labelled <- nrow(Y_labelled)
  n_unlabelled <- nrow(Y_unlabelled)
  n <- n_labelled + n_unlabelled
  p <- ncol(Y_labelled)
  
  # Initialization
  pi_k <- if (!is.null(pi_init)) pi_init else rep(1 / g, g)
  mu_k <- if (!is.null(mu_init)) mu_init else lapply(1:g, function(k) colMeans(Y_labelled[Z_labelled == k, , drop = FALSE]))
  
  # Determine covariance mode from Sigma_init
  cov_mode <- "shared"  # default
  if (!is.null(Sigma_init)) {
    if (is.matrix(Sigma_init)) {
      cov_mode <- "shared"
      Sigma_shared <- Sigma_init
      Sigma_k <- replicate(g, Sigma_shared, simplify = FALSE)
    } else if (is.list(Sigma_init) && length(Sigma_init) == g) {
      cov_mode <- "class"
      Sigma_k <- Sigma_init
      Sigma_shared <- NULL
    } else if (is.array(Sigma_init) && length(dim(Sigma_init)) == 3L && dim(Sigma_init)[3] == g) {
      cov_mode <- "class"
      Sigma_k <- lapply(1:g, function(i) Sigma_init[, , i, drop = FALSE][, , 1])
      Sigma_shared <- NULL
    } else {
      stop("`Sigma_init` must be a p times p matrix, a list of g p times p matrices, or a p times p times g array.")
    }
  } else {
    # ame covariance matrix estimated over all labeled data
    Sigma_shared <- cov(Y_labelled)
    Sigma_k <- replicate(g, Sigma_shared, simplify = FALSE)
  }
  
  alpha_k <- alpha_init
  
  # Combine data
  Y_all <- rbind(Y_labelled, Y_unlabelled)
  m_j <- c(rep(0, n_labelled), rep(1, n_unlabelled))
  Z_all <- c(Z_labelled, rep(NA, n_unlabelled))
  
  # Compute initial entropy
  Sigma_kk <- simplify2array(Sigma_k)
  current_ests <- list(pi = pi_k, mu = do.call(cbind, mu_k), sigma = Sigma_kk)
  d2_yj <- log(get_entropy(dat = Y_all, n = n, p = p, g = g, paralist = current_ests))
  
  # Logistic model initialization
  init_df <- data.frame(m_is_missing = m_j, log_entropy = d2_yj)
  glm_init <- glm(m_is_missing ~ log_entropy, data = init_df, family = binomial(link = "logit"))
  xi0 <- coef(glm_init)[1]
  xi1 <- coef(glm_init)[2]
  
  for (iter in 1:warm_up_iter) {
    alpha_old <- alpha_k 
    # =========================
    # E-Step
    # =========================
    gamma <- matrix(0, nrow = n, ncol = g)
    for (j in 1:n) {
      if (m_j[j] == 0) {
        k <- Z_all[j]
        gamma[j, k] <- 1
      } else {
        log_probs <- sapply(1:g, function(i) {
          if (cov_mode == "shared") {
            log(pi_k[i]) + dmvnorm(Y_all[j, ], mean = mu_k[[i]], sigma = Sigma_k[[i]], log = TRUE)
          } else {
            log(pi_k[i]) + dmvnorm(Y_all[j, ], mean = mu_k[[i]], sigma = Sigma_k[[i]], log = TRUE)
          }
        })
        gamma[j, ] <- normalise_logprob(log_probs)
      }
    }
    
    # =========================
    # M-Step (update parameters)
    # =========================
    Nk <- colSums(gamma)
    pi_k <- Nk / n
    
    for (i in 1:g) {
      mu_k[[i]] <- colSums(gamma[, i] * Y_all) / Nk[i]
    }
    
    if (cov_mode == "shared") {
      # Estimate shared Sigma (weighted over all components)
      Sigma_sum <- matrix(0, p, p)
      for (j in 1:n) {
        for (i in 1:g) {
          diff <- as.numeric(Y_all[j, ]) - mu_k[[i]]
          Sigma_sum <- Sigma_sum + gamma[j, i] * (diff %*% t(diff))
        }
      }
      Sigma_shared <- Sigma_sum / n
      Sigma_k <- replicate(g, Sigma_shared, simplify = FALSE)
    } else {
      # Class-specific Sigma_i
      for (i in 1:g) {
        Sigma_sum_i <- matrix(0, p, p)
        for (j in 1:n) {
          diff <- as.numeric(Y_all[j, ]) - mu_k[[i]]
          Sigma_sum_i <- Sigma_sum_i + gamma[j, i] * (diff %*% t(diff))
        }
        Sigma_k[[i]] <- Sigma_sum_i / Nk[i]
      }
    }
    
    # ========== entropy + m1j/m2j + xi updates ==========
    Sigma_kk <- simplify2array(Sigma_k)
    current_ests <- list(pi = pi_k, mu = do.call(cbind, mu_k), sigma = Sigma_kk)
    d2_yj <- log(get_entropy(dat = Y_all, n = n, p = p, g = g, paralist = current_ests))
    q_yj <- plogis(xi0 + xi1 * d2_yj)
    
    # Missing mechanism expectations
    m1j_k <- m2j_k <- numeric(n)
    for (j in 1:n) {
      if (m_j[j] == 1) {
        denom <- alpha_k + (1 - alpha_k) * q_yj[j]
        m1j_k[j] <- alpha_k / denom
        m2j_k[j] <- (1 - alpha_k) * q_yj[j] / denom
      }
    }
    
    # Update alpha
    alpha_k <- mean(m1j_k) 
    
    # Update xi via weighted logistic regression
    logit_df <- data.frame(log_entropy = d2_yj, m2j_k = m2j_k)
    glm_fit <- glm(m2j_k ~ log_entropy, data = logit_df, family = binomial(), weights = rep(1, nrow(logit_df)))
    xi0 <- coef(glm_fit)[1]; xi1 <- coef(glm_fit)[2]
    
    cat(sprintf("Iter %d:  alpha=%.4f | xi0=%.4f | xi1=%.4f | sum(m2j_k)=%.0f / %.0f\n",
                iter, alpha_k, xi0, xi1, sum(m2j_k), sum(m_j)))
    
    if (abs(alpha_k - alpha_old) < 1e-6) {
      cat("Early stopping: alpha converged.\n")
      break
    }
  }
  
  # Return Sigma in the same structure: shared matrix or list of g matrices
  Sigma_out <- if (cov_mode == "shared") Sigma_k[[1]] else Sigma_k
  
  return(list(pi = pi_k, mu = mu_k, Sigma = Sigma_out, xi = c(xi0, xi1),
              alpha = alpha_k, gamma = gamma, d2_yj = d2_yj)) 
}

