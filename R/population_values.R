#' @rdname fleiss
#' @export
bp_pop <- function(mu, sigma, values, kind = 1) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  mean_diff <- (mean(mu^2) - mean(mu)^2) * r
  cov_sum <- sum(sigma) / r
  c1 <- bp_aggr_get_c1(values, kind)
  d <- 2 / (r - 1) * (trace + mean_diff - cov_sum)
  1 - d / c1
}

#' @rdname fleiss
#' @export
conger_pop <- function(mu, sigma) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  top <- sum(sigma) - trace
  bottom <- (r - 1) * trace + r^2 * (mean(mu^2) - mean(mu)^2)
  top / bottom
}

#' @rdname fleiss
#' @export
cohen_pop <- function(mu, sigma) conger_pop(mu, sigma)

#' @rdname fleiss
#' @export
fleiss_pop <- function(mu, sigma) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  top <- sum(sigma) - trace - r * (mean(mu^2) - mean(mu)^2)
  bottom <- (r - 1) * trace + (r - 1) * r * (mean(mu^2) - mean(mu)^2)
  top / bottom
}
