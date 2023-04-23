#' Asymptotic variance for Brennan-Prediger.
#'
#' @param x Covariance matrix as calculated by `cov_mat_kappa`.
#' @param type Either `adf`, `normal`, or `elliptical`.
#' @param c1 Denominator in BP.
#' @return The asymptotic variance for the kappas.
#' @keywords internal
avar_bp <- \(x, type, c1) {
  p <- colSums(!is.na(x)) / nrow(x)
  mu <- colMeans(x, na.rm = TRUE)
  sigma <- stats::cov(x, use = "pairwise.complete.obs")
  gamma <- gamma_est(x, sigma, type)
  var <- cov_mat_bp(p, mu, sigma, gamma, x)
  avar_bp_(var, mu, sigma, c1)
}

cov_mat_bp <- function(p, mu, sigma, gamma, x) {
  r <- length(p)

  # The variance of mean(mu^2) - mean(mu)^2
  p_mat <- diag(1 / p)
  mu_middle <- (sigma + (p_mat - diag(r)) * diag(sigma))
  mu_vec <- (diag(r) - matrix(1, r, r) / r) %*% mu
  var_a <- c(4 * t(mu_vec) %*% mu_middle %*% mu_vec)

  # The covariances involving Sigma.
  gamma_pi <- gamma * pi_mat_empirical(x)
  d <- get_diag_indices(r, vech = TRUE)
  ones <- rep(1, length(d))
  d_minus_ones <- d*(1 + 1/r) - 2/r * ones
  var_b <- c(t(d_minus_ones) %*% gamma_pi %*% d_minus_ones)
  var_a + var_b
}

avar_bp_ <- function(vars, mu, sigma, c1) {
  r <- nrow(sigma)
  4 / (r - 1)^2 / c1^2 * vars
}
