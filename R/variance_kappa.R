#' Asymptotic variance for the kappas.
#'
#' @param x Covariance matrix as calculated by `cov_mat_kappa`.
#' @param type Either `adf`, `normal`, or `elliptical`.
#' @param fleiss If `TRUE` calculates variance for Fleiss kappa, else
#'    Conger's kappa.
#' @return The asymptotic variance for the kappas.
#' @keywords internal
avar <- \(x, type, fleiss) {
  p <- colSums(!is.na(x)) / nrow(x)
  mu <- colMeans(x, na.rm = TRUE)
  sigma <- stats::cov(x, use = "pairwise.complete.obs")
  gamma <- gamma_est(x, sigma, type)
  mat <- cov_mat_kappa(p, mu, sigma, gamma, x)
  avar_(mat, mu, sigma, fleiss)
}

#' Asymptotic variance for the kappas.
#'
#' @param mat Covariance matrix as calculated by `cov_mat_kappa`.
#' @param mu Vector of means.
#' @param sigma Covariance matrix
#' @param fleiss If `TRUE` calculates variance for Fleiss kappa, else
#'    Conger's kappa.
#' @return The asymptotic variance for the kappas.
#' @keywords internal
avar_ <- function(mat, mu, sigma, fleiss) {
  r <- nrow(sigma)
  yy <- sum(diag(sigma))
  xx <- sum(sigma) - yy
  zz <- mean(mu^2) - mean(mu)^2

  vec <- if (fleiss) {
    mult <- 1 / ((r - 1) * (r * zz + yy)^2)
    c(r * zz + yy, r * zz - xx, -r * (xx + yy))
  } else {
    mult <- 1 / (r^2 * zz + (r - 1) * yy)^2
    c(r^2 * zz + (r - 1) * yy, -(r - 1) * xx, -r^2 * xx)
  }

  c(t(vec) %*% mat %*% vec) * mult^2
}

#' Covariance matrix (x,y,z) for Fleiss' and Cohen's kappa.
#'
#' Returns the covariance matrix of x, y, z as defined in the paper on
#'   quadratic Fleiss kappa and Cohen's kappa.
#'
#' @param p Vector of probabilities for being missing.
#' @param mu Vector of means.
#' @param sigma Covariance matrix
#' @param gamma Covariance matrix of the covariances.
#' @return The covariance matrix of x, y, z
#' @keywords internal
cov_mat_kappa <- function(p, mu, sigma, gamma, x) {
  r <- length(p)

  # The covariances involving s.
  gamma_pi <- gamma * pi_mat_empirical(x)
  d <- get_diag_indices(r, vech = TRUE)
  ones <- rep(1, length(d))
  ones_minus_d <- ones - d
  cov_ss_ss <- t(ones_minus_d) %*% (gamma_pi) %*% ones_minus_d
  cov_ss_tr <- t(ones_minus_d) %*% (gamma_pi) %*% d
  cov_tr_tr <- t(d) %*% gamma_pi %*% d

  # The variance of mean(mu^2) - mean(mu)^2
  p_mat <- diag(1 / p)
  mu_middle <- (sigma + (p_mat - diag(r)) * diag(sigma))
  mu_vec <- (diag(r) - matrix(1 / r, r, r)) %*% mu
  cov_mu_mu <- c(4 / r^2 * t(mu_vec) %*% mu_middle %*% mu_vec)
  matrix(c(
    4 * cov_ss_ss, 2 * cov_ss_tr, 0,
    2 * cov_ss_tr, cov_tr_tr, 0,
    0, 0, cov_mu_mu
  ), nrow = 3, byrow = FALSE)
}
