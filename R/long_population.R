bp_pop <- function(mu, sigma, c1) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  mean_diff <- (mean(mu^2) - mean(mu)^2) * r
  cov_sum <- sum(sigma) / r
  d <- 2 / (r - 1) * (trace + mean_diff - cov_sum)
  1 - d / c1
}

conger_pop <- function(mu, sigma) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  top <- sum(sigma) - trace
  bottom <- (r - 1) * trace + r^2 * (mean(mu^2) - mean(mu)^2)
  top / bottom
}

cohen_pop <- function(mu, sigma) conger_pop(mu, sigma)

fleiss_pop <- function(mu, sigma) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  top <- sum(sigma) - trace - r * (mean(mu^2) - mean(mu)^2)
  bottom <- (r - 1) * trace + (r - 1) * r * (mean(mu^2) - mean(mu)^2)
  top / bottom
}
