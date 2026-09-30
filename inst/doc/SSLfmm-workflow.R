## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----simulation---------------------------------------------------------------
library(SSLfmm)

mu <- matrix(c(-1.5, 1.5), nrow = 1, ncol = 2)
sim <- simulate_mixed_missingness(
  n = 120,
  pi = c(0.5, 0.5),
  mu = mu,
  sigma = matrix(1, 1, 1),
  alpha = 0.10,
  mar_rate = 0.25,
  seed = 2026
)

head(sim$data)
table(sim$data$missing_source, useNA = "ifany")

## ----fit----------------------------------------------------------------------
x <- as.matrix(sim$data["x1"])
fit <- fit_sslfmm(
  x,
  sim$data$label,
  g = 2,
  method = "mixed",
  covariance_type = "equal",
  indicator = "unknown",
  n_starts = 3,
  seed = 2027
)

fit
summary(fit)

## ----prediction---------------------------------------------------------------
pred_class <- predict(fit, x)
posterior <- predict(fit, x, type = "posterior")
entropy <- predict(fit, x, type = "entropy")

head(pred_class)
head(posterior)
head(entropy)

## ----assessment---------------------------------------------------------------
perf <- classification_performance(
  sim$data$truth,
  pred_class,
  posterior
)

perf$metrics
perf$confusion_matrix

## ----data---------------------------------------------------------------------
data("blood_transfusion")
head(blood_transfusion)
table(blood_transfusion$missing_indicator)

