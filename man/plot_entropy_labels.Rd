\name{plot_entropy_labels}
\alias{plot_entropy_labels}
\title{Boxplot of Posterior Entropy by Class Label}
\usage{
plot_entropy_labels(entropy, labels = NULL,
  main = "Posterior entropy by class label",
  xlab = "Class label", ylab = "Posterior entropy",
  legend = TRUE, ...)
}
\arguments{
\item{entropy}{Numeric posterior-entropy vector or an \code{SSLfmm} fit.}
\item{labels}{Class labels. If an \code{SSLfmm} fit is supplied, fitted classes are used by default.}
\item{main,xlab,ylab}{Plot labels.}
\item{legend}{Retained for backward compatibility; boxplots do not require a legend.}
\item{...}{Additional arguments passed to \code{graphics::boxplot()}.}
}
\value{Invisibly returns the legacy entropy-group summaries together with the boxplot summary.}
\description{
Draws a boxplot of finite posterior-entropy values grouped by class label.
}
\examples{
ent <- c(0.05, 0.10, 0.20, 0.40, 0.55, 0.65)
lab <- c(1, 1, 2, 1, 2, 2)
plot_entropy_labels(ent, lab)
}
