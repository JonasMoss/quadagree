### ============================================================================
###   Testing inference for Brennan-Prediger.
###
###  1. Make an estimator function.
###  2. Deduce and implement the variance.
###  3. Refactor to make speedier.
###  4. Implement simple solution in the package.
###  5. Make it a "fleissci" element.
### ============================================================================

x <- dat.zapf2016

bp2 <- \(x) {
  y <- as.matrix(x)
}

bp_pop <- \(sigma, mu, type = 1) {
  stopifnot(type == 1 || type == 2)
  if(type == 1) {

  }

  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  mean_diff <- (mean(mu^2) - mean(mu)^2)
  top <- sum(sigma) + r * mean_diff - sum(sigma) / r
  1 - 2 / (r - 1) / c1
}
