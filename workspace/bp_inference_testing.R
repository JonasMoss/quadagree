### ============================================================================
###
###  Testing the inference formulas for BP.
###
### ============================================================================

bp_est <- \(x, values = seq(min(x), max(x))) {
  n <- nrow(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs") * (n - 1) / n
  mu <- colMeans(x, na.rm = TRUE)
  bp_pop(mu, sigma, values)
}

bp_pop <- \(mu, sigma, values, type = 1) {
  r <- ncol(sigma)

  c1 <- bp_get_c1(values, type)
  trace <- sum(diag(sigma))
  top <- sum(sigma) - trace / r + (mean(mu^2) - mean(mu)^2)

  bottom <-
  top / bottom
}
