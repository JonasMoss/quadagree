#' Calculate the quadratically weighted Fleiss's kappa, Conger's kappa, and
#'    Cohen's kappa
#'
#' @param x Input data on wide. Columns are judges, items are ratings. Missing
#'    values are allowed.
#' @param mu Population vector of means.
#' @param sigma Population covariance matrix.
#' @return Sample quadratically weighted Fleiss' kappa.
#' @examples
#' x <- irrCAC::cac.raw4raters[2:9, ]
#' fleiss(x)
#' # [1] 0.6666667
#' irrCAC::fleiss.kappa.raw(x, weights = "quadratic")$est$coeff.val
#' # [1] 0.66667 (irrCAC prematurely rounds the coefficient)
#'
#' x <- irrCAC::cac.raw4raters[2:9, ]
#' conger(x)
#' # [1] 0.6719243
#' irrCAC::conger.kappa.raw(x, weights = "quadratic")$est$coeff.val
#' # [1] 0.67192 (irrCAC prematurely rounds the coefficient)
#' @export

#' @name fleiss
fleiss <- function(x) {
  n <- nrow(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs") * (n - 1) / n
  mu <- colMeans(x, na.rm = TRUE)
  fleiss_pop(mu, sigma)
}

#' @rdname fleiss
#' @export
conger <- function(x) {
  n <- nrow(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs") * (n - 1) / n
  mu <- colMeans(x, na.rm = TRUE)
  conger_pop(mu, sigma)
}

#' @rdname fleiss
#' @export
cohen <- function(x) conger_pop(x)

#' Calculates the population value of Conger's kappa and Fleiss' kappa.
#'
#' @param mu Vector of means.
#' @param sigma Covariance matrix; must have the same dimensions as mu.
#' @return Population value of Conger's kappa / Fleiss' kappa.

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
