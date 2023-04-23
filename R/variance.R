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

#' Get the indices of the diagonal for the variances.
#' @param r Number of raters.
#' @param vech If `TRUE`, returns for the result for half-vectorization.
#' @return A vector of `TRUE` where the variances are.
#' @keywords internal
get_diag_indices <- function(r, vech = TRUE) {
  e_mat <- matrixcalc::elimination.matrix(r)
  indices <- rep(0, r^2)
  indices[c(1, (r + 1) * (1:(r - 1)) + 1)] <- 1
  if (vech) c(e_mat %*% indices) else c(indices)
}

#' Gamma matrix
#'
#' Calculate the gamma matrix from a matrix of observations.
#' @param x A numeric matrix of observations.
#' @param sigma Covariance matrix of the data.
#' @param type One of `adf`, `normal` and `elliptical`.
#' @return The sample estimate of the gamma matrix.
#' @keywords internal
gamma_est <- function(x, sigma, type = "adf") {
  if (type == "adf") {
    i_row <- \(n) unlist(lapply(seq_len(n), seq.int, n))
    i_col <- \(n) rep.int(seq_len(n), times = rev(seq_len(n)))
    rows <- i_row(ncol(x))
    cols <- i_col(ncol(x))
    y <- t(x) - colMeans(x, na.rm = TRUE)
    z <- y[cols, , drop = FALSE] * y[rows, , drop = FALSE]
    mat <- z - rowMeans(z, na.rm = TRUE)
    nas <- is.na(mat)
    mat[nas] <- 0
    return(base::tcrossprod(mat) / base::tcrossprod(!nas))
  }

  k <- ncol(sigma)
  e_mat <- matrixcalc::elimination.matrix(k)
  k_mat <- matrixcalc::K.matrix(k)
  gamma <- ((diag(k^2) + k_mat) %*% (sigma %x% sigma))
  gamma <- gamma * kurtosis_correction(x, type = type)
  e_mat %*% gamma %*% t(e_mat)
}

#' Calculate unbiased sample kurtosis.
#' @param x Matrix of valus.
#' @return Unbiased sample kurtosis.
#' @keywords internal
kurtosis <- function(x) {
  n <- nrow(x)
  g2 <- \(x) mean((x - mean(x, na.rm = TRUE))^4) / stats::var(x, na.rm = TRUE)^2
  kurtosis <- \(x) (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2(x) + 6)
  mean(apply(x, 2, kurtosis), na.rm = TRUE) - 3
}

#' Calculate kurtosis correction
#' @param x Matrix of values
#' @param type The type of correction, either "normal" or "elliptical".
#' @keywords internal
kurtosis_correction <- function(x, type) {
  kurt <- if (type == "normal") 0 else kurtosis(x)
  1 + kurt / 3
}

#' Calculate the capital pi matrix.
#'
#' The capital pi matrix is multiplied element-wise with the asymptotic
#'    covariance matrix of s to correct for missing values.
#'
#' @param p Vector of probabilities for being missing.
#' @param vech Does not affect anything. The half-vectorized matrix is always
#'    returned.
#' @return The capital pi matrix.
#' @keywords internal
pi_mat <- function(p, vech = TRUE) {
  r <- length(p)
  f <- \(x) prod(p[unique(x)])
  indices <- arrangements::combinations(r, 2, replace = TRUE)
  ps <- apply(indices, 1, f)
  g <- Vectorize(\(x, y) f(c(indices[x, ], indices[y, ])) / (ps[x] * ps[y]))
  x <- seq_len(nrow(indices))
  outer(x, x, g)
}

#' Calculates the empirical capital pi matrix.
#'
#' @param x Vector of probabilities for being missing.
#' @return The capital pi matrix.
#' @keywords internal
pi_mat_empirical <- \(x) {
  r <- ncol(x)
  ind2 <- arrangements::combinations(r, 2, replace = TRUE)
  ind4 <- arrangements::combinations(seq_len(nrow(ind2)), 2, replace = TRUE)

  p2_hats <- apply(ind2, 1, \(i) mean(!is.na(x[, i[1]]) & !is.na(x[, i[2]])))
  hats <- apply(ind4, 1, \(i) {
    j1 <- ind2[i[1], ]
    j2 <- ind2[i[2], ]
    mean(!is.na(x[, j1[1]]) & !is.na(x[, j1[2]]) &
      !is.na(x[, j2[1]]) & !is.na(x[, j2[2]])) /
      (p2_hats[i[1]] * p2_hats[i[2]])
  })

  new_mat <- matrix(NA, choose(r + 1, 2), choose(r + 1, 2))
  new_mat[ind4] <- hats
  as.matrix(Matrix::forceSymmetric(new_mat, uplo = "U"))
}
