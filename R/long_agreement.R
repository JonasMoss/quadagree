fleiss <- function(x) {
  n <- nrow(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(x, na.rm = TRUE)
  fleiss_pop(mu, sigma)
}

bp <- function(x, c1) {
  y <- as.matrix(x)
  n <- nrow(y)
  sigma <- stats::cov(y, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(y, na.rm = TRUE)
  bp_pop(mu, sigma, c1)
}

conger <- function(x) {
  n <- nrow(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) stop("The data does not contain sufficient non-NAs.")
  mu <- colMeans(x, na.rm = TRUE)
  conger_pop(mu, sigma)
}

cohen <- function(x) conger(x)

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
