#' Normalise Log-Probabilities
#'
#' Converts log-probabilities into a probability distribution
#' by exponentiating in a numerically stable way.
#'
#' @param log_probs A numeric vector of log-probabilities.
#' @return A numeric vector summing to 1 (the normalised probabilities).
#' @examples
#' lp <- c(-1000, -999, -998)
#' normalise_logprob(lp)
#' @export
normalise_logprob <- function(log_probs) {
  stabilized_probs <- exp(log_probs - logsumexp(log_probs))
  return(stabilized_probs)
}