.safe_div <- function(a, b, zero = NA_real_) {
  ifelse(b > 0, a / b, zero)
}

#' Classification performance metrics
#'
#' @param truth True class labels.
#' @param predicted Predicted class labels.
#' @param posterior Optional matrix of class posterior probabilities, one row per observation.
#' @return A list with overall metrics, per-class metrics and a confusion matrix.
#' @details Macro averages are taken over classes represented in `truth` after
#'   removing incomplete truth/prediction pairs. If a represented true class is
#'   never predicted, its precision and F1 score are defined as zero rather than
#'   omitted from the macro average.
#' @export
classification_performance <- function(truth, predicted, posterior = NULL) {
  if (length(truth) != length(predicted)) stop("'truth' and 'predicted' must have the same length.", call. = FALSE)
  keep <- !is.na(truth) & !is.na(predicted)
  if (!any(keep)) stop("No complete truth/prediction pairs are available.", call. = FALSE)
  tt <- as.character(truth[keep])
  pp <- as.character(predicted[keep])
  lev <- sort(unique(c(tt, pp)))
  cm <- table(factor(tt, levels = lev), factor(pp, levels = lev),
              dnn = c("truth", "predicted"))
  tp <- diag(cm)
  support <- rowSums(cm)
  predicted_n <- colSums(cm)

  # Macro metrics are defined over classes that occur in the truth labels.
  # A true class that receives no predictions must contribute zero precision
  # and zero F1, rather than being silently removed with na.rm = TRUE.
  truth_class <- support > 0
  recall <- .safe_div(tp, support)
  precision <- .safe_div(tp, predicted_n, zero = 0)
  f1 <- .safe_div(2 * precision * recall, precision + recall, zero = 0)

  accuracy <- sum(tp) / sum(cm)
  metrics <- c(accuracy = accuracy,
               error_rate = 1 - accuracy,
               balanced_accuracy = mean(recall[truth_class]),
               macro_precision = mean(precision[truth_class]),
               macro_recall = mean(recall[truth_class]),
               macro_f1 = mean(f1[truth_class]))

  if (!is.null(posterior)) {
    posterior <- as.matrix(posterior)
    if (nrow(posterior) != length(truth)) stop("'posterior' must have one row per observation.", call. = FALSE)
    post <- posterior[keep, , drop = FALSE]
    if (!is.numeric(post) || any(!is.finite(post))) {
      stop("'posterior' must contain only finite numeric probabilities.", call. = FALSE)
    }
    if (any(post < 0)) {
      stop("'posterior' probabilities must be nonnegative.", call. = FALSE)
    }
    if (is.null(colnames(post))) {
      if (ncol(post) != length(lev)) {
        stop("Without column names, 'posterior' must have one column per class represented by 'truth'/'predicted'.", call. = FALSE)
      }
      colnames(post) <- lev
    } else if (anyDuplicated(colnames(post))) {
      stop("'posterior' column names must be unique class labels.", call. = FALSE)
    }
    idx <- match(tt, colnames(post))
    if (anyNA(idx)) stop("'posterior' must contain a column for every true class.", call. = FALSE)
    rs <- rowSums(post)
    if (any(!is.finite(rs)) || any(rs <= 0)) stop("Each 'posterior' row must contain positive finite probability mass.", call. = FALSE)
    post <- post / rs
    ptrue <- post[cbind(seq_along(idx), idx)]
    log_loss <- -mean(log(pmax(ptrue, 1e-15)))
    onehot <- matrix(0, nrow(post), ncol(post))
    onehot[cbind(seq_along(idx), idx)] <- 1
    brier <- mean(rowSums((post - onehot)^2))
    metrics <- c(metrics, log_loss = log_loss, brier_score = brier)
  }

  per_class <- data.frame(class = lev, support = as.integer(support),
                          precision = as.numeric(precision), recall = as.numeric(recall),
                          f1 = as.numeric(f1), row.names = NULL)
  list(metrics = metrics, per_class = per_class, confusion_matrix = cm)
}

.entropy_groups <- function(entropy) {
  n <- length(entropy)
  ord <- order(entropy, na.last = NA)
  grp <- rep(NA_character_, n)
  if (!length(ord)) return(factor(grp, levels = c("Low", "Medium", "High")))
  rank_index <- seq_along(ord)
  cut1 <- ceiling(length(ord) / 3)
  cut2 <- ceiling(2 * length(ord) / 3)
  lab <- ifelse(rank_index <= cut1, "Low", ifelse(rank_index <= cut2, "Medium", "High"))
  grp[ord] <- lab
  factor(grp, levels = c("Low", "Medium", "High"))
}

#' Boxplot of posterior entropy by class label
#'
#' Displays the distribution of posterior entropy separately for each class label. If an
#' SSLfmm fit is supplied, fitted entropy and predicted classes are used by default.
#'
#' @param entropy Numeric entropy vector or an SSLfmm object.
#' @param labels Optional class labels. With an SSLfmm object, predicted classes are used when omitted.
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param legend Retained for backward compatibility; boxplots do not require a legend.
#' @param ... Additional arguments passed to graphics::boxplot().
#' @return Invisibly, a list containing entropy_group, counts, proportions, and the boxplot summary.
#' @export
plot_entropy_labels <- function(entropy, labels = NULL,
                                main = "Posterior entropy by class label",
                                xlab = "Class label", ylab = "Posterior entropy",
                                legend = TRUE, ...) {
  if (inherits(entropy, "SSLfmm")) {
    object <- entropy
    entropy <- object$entropy
    if (is.null(labels)) {
      idx <- max.col(object$posterior, ties.method = "first")
      labels <- object$label_levels[idx]
    }
  }
  entropy <- as.numeric(entropy)
  if (is.null(labels)) stop("'labels' is required when 'entropy' is a numeric vector.", call. = FALSE)
  if (length(labels) != length(entropy)) stop("'labels' and 'entropy' must have the same length.", call. = FALSE)
  ok <- is.finite(entropy) & !is.na(labels)
  if (!any(ok)) stop("No finite entropy values with non-missing labels are available.", call. = FALSE)

  labs <- factor(as.character(labels))
  grp <- .entropy_groups(entropy)
  counts <- table(labs[ok], grp[ok])
  proportions <- sweep(counts, 2L, pmax(colSums(counts), 1), "/")

  bp <- graphics::boxplot(entropy[ok] ~ droplevels(labs[ok]),
                          main = main, xlab = xlab, ylab = ylab, ...)
  invisible(list(entropy_group = grp, counts = counts, proportions = proportions,
                 boxplot = bp))
}
