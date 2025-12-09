#' Numerically Stable Log-Sum-Exp
#'
#' Computes \eqn{\log(\sum_i \exp(x_i))} in a numerically stable way
#' by subtracting the maximum value before exponentiation.
#'
#' @param x A numeric vector.
#'
#' @return A single numeric value: the log-sum-exp of `x`.
#' @examples
#' logsumexp(c(1000, 1001, 1002))
#' @export
logsumexp <- function(x) {
  max_x <- max(x)
  stabilized_exp <- exp(x - max_x)
  log_sum_exp <- log(sum(stabilized_exp)) + max_x
  return(log_sum_exp)
}
