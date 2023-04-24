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
