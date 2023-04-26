#' Calculate limits of a confidence interval.
#'
#' @param alternative Alternative choosen.
#' @param conf_level Confidence level.
#' @keywords internal
limits <- function(alternative, conf_level) {
  half <- (1 - conf_level) / 2
  if (alternative == "two.sided") {
    return(c(half, 1 - half))
  }
  if (alternative == "greater") {
    return(c(2 * half, 1))
  }
  if (alternative == "less") {
    return(c(0, conf_level))
  }
}

#' Get the indices of the diagonal for the variances.
#' @param r Number of raters.
#' @param vech If `TRUE`, returns for the result for half-vectorization.
#' @return A vector of `TRUE` where the variances are.
#' @keywords internal
get_diag_indices <- function(r, vech = TRUE) {
  if (vech) {
    indices <- rep(0, choose(r + 1, 2))
    indices[c(1, 1 + cumsum(r:2))] <- 1
  } else {
    indices <- rep(0, r^2)
    indices[c(1, (r + 1) * (1:(r - 1)) + 1)] <- 1
  }
  indices
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
    z <- y[cols, ] * y[rows, ]
    mat <- z - rowMeans(z, na.rm = TRUE)
    if (!anyNA(mat)) {
      div <- nrow(x)
    } else {
      nas <- is.na(mat)
      mat[nas] <- 0
      div <- tcrossprod(!nas)
    }
    return(tcrossprod(mat) / div)
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
  if (!anyNA(x)) {
    return(matrix(1, choose(r + 1, 2), choose(r + 1, 2)))
  }

  ind2 <- arrangements::combinations(r, 2, replace = TRUE)
  ind4 <- arrangements::combinations(seq_len(nrow(ind2)), 2, replace = TRUE)

  nisna <- !is.na(x)
  combs <- apply(ind2, 1, \(i) nisna[, i[1]] & nisna[, i[2]])
  p2_hats <- colMeans(combs)

  hats <- apply(ind4, 1, \(i) {
    mean(combs[, i[1]] & combs[, i[2]]) / (p2_hats[i[1]] * p2_hats[i[2]])
  })

  new_mat <- matrix(NA, choose(r + 1, 2), choose(r + 1, 2))
  new_mat[ind4] <- hats
  as.matrix(Matrix::forceSymmetric(new_mat, uplo = "U"))
}

#' Get the required value of c1.
#' @param values Vector of values.
#' @param kind The kind of c1 requested.
#' @return The value of c1.
#' @keywords internal
bp_get_c1 <- \(values, kind) {
  if (kind == 1) {
    n_cat <- length(values)
    (2 * n_cat * sum(values^2) - 2 * sum(values)^2) / n_cat^2
  } else {
    0.5 * (max(values) - min(values))^2
  }
}

#' @export
print.quadagree <- function(x, digits = getOption("digits"), ...) {
  at <- \(y) attr(x, y)
  cat("Call: ", paste(deparse(at("call")),
    sep = "\n",
    collapse = "\n"
  ), "\n\n", sep = "")

  if (!is.null(x)) {
    cat(format(100 * at("conf_level")),
      "% confidence interval (n = ", at("n"), ").\n",
      sep = ""
    )
    print(x[1:2], digits = digits)
    cat("\n")
  }

  if (!is.null(at("estimate"))) {
    cat("Sample estimates.\n")
    print(
      c(
        kappa = at("estimate"),
        sd = at("sd")
      ),
      digits = digits
    )
  }
  invisible(x)
}

get_transformer <- function(transform) {
  transformers <- list(
    fisher = transformer_fisher,
    none = transformer_none,
    arcsin = transformer_arcsin,
    log = transformer_log
  )

  if (transform %in% names(transformers)) {
    transformers[[transform]]
  } else {
    stop(paste0("`transformer = ", transform, "` not supported."))
  }
}

transformer_fisher <- c(
  est = \(est) atanh(est),
  sd = \(est, sd) sd / (1 - est^2),
  inv = tanh
)

transformer_log <- c(
  est = \(est) log(1 - est),
  sd = \(est, sd) sd / abs(1 - est),
  inv = \(x) 1 - exp(x)
)

transformer_none <- c(
  est = \(est) est,
  sd = \(est, sd) sd,
  inv = \(x) x
)

transformer_arcsin <- c(
  est = asin,
  sd = \(est, sd) sd / sqrt(1 - est^2),
  inv = sin
)

tr <- \(x) {
  sum(x[1L + 0L:(dim(x)[1L] - 1L) * (dim(x)[1L] + 1L)])
}
