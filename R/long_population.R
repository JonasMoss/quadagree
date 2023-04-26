bp_pop <- function(mu, sigma, c1) {
  r <- ncol(sigma)
  trace <- tr(sigma)
  mean_diff <- (mean(mu^2) - mean(mu)^2) * r
  cov_sum <- sum(sigma) / r
  d <- 2 / (r - 1) * (trace + mean_diff - cov_sum)
  1 - d / c1
}

conger_pop <- function(mu, sigma) {
  r <- ncol(sigma)
  trace <- tr(sigma)
  mean_diff <- (mean(mu^2) - mean(mu)^2) * r^2
  top <- sum(sigma) - trace
  bottom <- (r - 1) * trace + mean_diff
  top / bottom
}

cohen_pop <- function(mu, sigma) conger_pop(mu, sigma)

fleiss_pop <- function(mu, sigma) {
  n <- length(mu)
  r <- ncol(sigma)
  trace <- tr(sigma)
  mean_diff <- (sum(mu^2) / n - (sum(mu) / n)^2) * r
  top <- sum(sigma) - trace - mean_diff
  bottom <- (r - 1) * (trace + mean_diff)
  top / bottom
}
