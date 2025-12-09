#' Complete-Data Warm-Up Initialization for Semi-Supervised FMM
#'   with a Mixed-Missingness Mechanism
#'
#' Uses both labeled and unlabeled subsets of the data to obtain quick initial
#' estimates for mixture parameters and missingness mechanism parameters
#' (\code{alpha}, \code{xi}) via a warm-up EM procedure.
#'
#' @param data A data frame containing:
#'   \itemize{
#'     \item The first \code{p} columns: numeric variables used in the FMM.
#'     \item A column \code{missing}: indicator (0 = labeled, 1 = unlabeled/missing).
#'     \item A column \code{z}: class labels for labeled rows (\code{1:g});
#'           \code{NA} for unlabeled.
#'   }
#' @param g Integer, number of mixture components (default \code{2}).
#' @param ncov Integer, covariance structure:
#'   \code{1} = shared (equal), \code{2} = class-specific (unequal).
#' @param alpha_init Numeric in (0,1), initial MCAR proportion (default \code{0.01}).
#' @param warm_up_iter Integer, number of warm-up EM iterations (default \code{200}).
#' @param tol Convergence tolerance on \code{alpha} (default \code{1e-6}).
#'
#' @return A list with initial values for \code{\link{EM_FMM_SemiSupervised}}:
#' \itemize{
#'   \item \code{pi}   - mixture weights.
#'   \item \code{mu}   - list of component mean vectors.
#'   \item \code{Sigma} - covariance matrix/matrices.
#'   \item \code{alpha} - MCAR proportion.
#'   \item \code{xi}    - logistic regression coefficients for MAR mechanism.
#' }
#'
#' @details
#' \itemize{
#'   \item This function first calls \code{\link{initialestimate}} to get initial
#'         \eqn{\pi}, \eqn{\mu}, \eqn{\Sigma}.
#'   \item Then it calls \code{\link{EM_FMM_SemiSupervised_Initial}} with these
#'         values for a short warm-up run.
#'   \item Covariance structure (\code{equal} vs. \code{unequal}) is determined
#'         by \code{ncov}.
#' }
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom matrixStats rowLogSumExps
#' @export
EM_FMM_SemiSupervised_Complete_Initial <- function(data, g = 2, ncov = 1,
                                                   alpha_init   = 0.01,
                                                   warm_up_iter = 200,
                                                   tol          = 1e-6) {
  
  data <- as.data.frame(data)
  p <- ncol(data) - 2L  # assumes last 2 cols are 'missing' and 'z'
  
  Y_labelled   <- data[data$missing == 0, seq_len(p), drop = FALSE]
  Y_unlabelled <- data[data$missing != 0, seq_len(p), drop = FALSE]
  Z_labelled   <- data[data$missing == 0, "z"]
  
  # Initial values from labeled subset
  inits <- initialestimate(dat = data[, seq_len(p), drop = FALSE],
                           zm  = data$z, g = g, ncov = ncov)
  pi_init    <- inits$pi
  mu_init    <- lapply(seq_len(ncol(inits$mu)), function(i) inits$mu[, i])
  Sigma_init <- inits$sigma
  
  # Warm-up EM to refine alpha and xi
  init_res <- EM_FMM_SemiSupervised_Initial(
    Y_labelled, Z_labelled, Y_unlabelled,
    g           = g,
    pi_init     = pi_init,
    mu_init     = mu_init,
    Sigma_init  = Sigma_init,
    alpha_init  = alpha_init,
    warm_up_iter = warm_up_iter,
    tol         = tol
  )
  
  list(
    pi    = init_res$pi,
    mu    = init_res$mu,
    Sigma = init_res$Sigma,
    alpha = init_res$alpha,
    xi    = init_res$xi
  )
}
