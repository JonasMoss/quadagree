avar <- \(x, type, fleiss) {
  p <- colSums(!is.na(x)) / nrow(x)
  mu <- colMeans(x, na.rm = TRUE)
  sigma <- stats::cov(x, use = "pairwise.complete.obs")
  gamma <- gamma_est(x, sigma, type)
  mat <- cov_mat(p, mu, sigma, gamma)
  avar_(mat, mu, sigma, fleiss)
}

#' Asymptotic variance for the kappas.
#'
#' @param mat Covariance matrix as calculated by `cov_mat`.
#' @param mu Vector of means.
#' @param sigma Covariance matrix
#' @param fleiss If `TRUE` calculates variance for Fleiss kappa, else
#'    Conger's kappa.
#' @return The asymptotic variance for the kappas.
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

#' The covariance matrix of the elements involved in the inference.
#'
#' Returns the covariance matrix of x, y, z as defined in the paper.
#'
#' @param p Vector of probabilities for being missing.
#' @param mu Vector of means.
#' @param sigma Covariance matrix
#' @param gamma Covariance matrix of the covariances.
#' @return The covariance matrix of x, y, z
#' @keywords internal
cov_mat <- function(p, mu, sigma, gamma) {
  r <- length(p)

  # The covariances involving s.
  gamma_pi <- gamma * pi_mat(p, vech = TRUE)
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
    4*cov_ss_ss, 2*cov_ss_tr, 0,
    2*cov_ss_tr, cov_tr_tr, 0,
    0, 0, cov_mu_mu
  ), nrow = 3, byrow = FALSE)
}

#' Get the indices of the diagonal for the variances.
#' @param r Number of raters.
#' @param vech If `TRUE`, returns for the result for half-vectorization.
#' @return A vector of `TRUE` where the variances are.
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
    n <- nrow(x)
    i_row <- \(n) unlist(lapply(seq_len(n), seq.int, n))
    i_col <- \(n) rep.int(seq_len(n), times = rev(seq_len(n)))
    rows <- i_row(ncol(x))
    cols <- i_col(ncol(x))
    y <- t(x) - colMeans(x, na.rm = TRUE)
    z <- y[cols, , drop = FALSE] * y[rows, , drop = FALSE]
    mat <- z - rowMeans(z, na.rm = TRUE)
    nas <- is.na(mat)
    mat[nas] <- 0
    base::tcrossprod(mat) / base::tcrossprod(!nas)
  } else {
    k <- ncol(sigma)
    e_mat <- matrixcalc::elimination.matrix(k)
    k_mat <- matrixcalc::K.matrix(k)
    gamma <- ((diag(k^2) + k_mat) %*% (sigma %x% sigma))
    gamma <- gamma * kurtosis_correction(x, type = type)
    e_mat %*% gamma %*% t(e_mat)
  }
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
#' @param vech If `TRUE`, returns for the result for half-vectorization.
#' @return The capital pi matrix.
pi_mat <- function(p, vech = TRUE) {
  r <- length(p)
  if (vech) {
    m <- (r * (r + 1)) / 2
    temp_mat <- matrix(rep(seq(r), r), nrow = r)
    ri <- c(matrixcalc::vech(temp_mat))
    ci <- c(matrixcalc::vech(t(temp_mat)))
  } else {
    m <- r^2
    temp_mat <- matrix(rep(1:r, r), nrow = r)
    ri <- c(matrixcalc::vec(temp_mat))
    ci <- c(matrixcalc::vec(t(temp_mat)))
  }
  outer(seq(m), seq(m), Vectorize(function(i, j) {
    numbers <- c()
    for (k in c(i, j)) {
      numbers <- c(numbers, ri[k])
      if (ri[k] != ci[k]) numbers <- c(numbers, ci[k])
    }
    table <- Rfast::Table(numbers)
    1 / prod(p[as.numeric(names(table))[table >= 2]])
  }))
}
