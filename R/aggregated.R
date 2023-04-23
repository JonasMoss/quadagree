bp_aggr_prepare <- \(x, values, type) {
  r <- sum(x[1, ])
  y <- as.matrix(x)
  calc <- bp_aggr_calc(y, values)
  c1 <- bp_aggr_get_c1(values, type)
  list(calc = calc, c1 = c1, r = r)
}

bp_aggr_calc <- \(y, values) {
  xtx <- c(tcrossprod(values^2, y))
  xt12 <- c(tcrossprod(values, y))^2
  cbind(xtx, xt12)
}

bp_aggr_est <- \(x, values, type) {
  args <- bp_aggr_prepare(x, values, type)
  do.call(bp_aggr_est_matrix, args)
}

bp_aggr_est_matrix <- \(calc, c1, r) {
  means <- colMeans(calc)
  disagreement <- 2 / (r - 1) * (means[1] - 1 / r * means[2])
  unname(1 - disagreement / c1)
}

bp_aggr_get_c1 <- \(values, type) {
  if (type == 1) {
    w <- outer(values, values, Vectorize(\(x, y) (x - y)^2))
    n_cat <- length(values)
    sum(w) / n_cat^2
  } else {
    0.5 * (max(values) - min(values))^2
  }
}

bp_aggr_var <- \(x, values, type) {
  args <- bp_aggr_prepare(x, values, type)
  do.call(bp_aggr_var_matrix, args)
}

bp_aggr_var_matrix <- \(calc, c1, r) {
  n <- nrow(calc)
  phi <- stats::cov(calc) * (n - 1) / n
  phi[1, 1] <- phi[1, 1]
  phi[1, 2] <- phi[2, 1] <- -phi[1, 2] / r
  phi[2, 2] <- phi[2, 2] / r^2
  1 / c1^2 * 4 / (r - 1)^2 * sum(phi)
}

#' Variance for Fleiss' kappa with aggregated raters.
#' @param x Data on Fleiss form.
#' @param values Values to attach to each column on the Fleiss form data.
#'    Defaults to `1:C`, where `C` is the number of categories.
#' @return Calculated value of Fleiss' kappa.
#' @keywords internal
fleiss_aggr_var <- \(x, values = seq_len(ncol(x))) {
  r <- sum(x[1, ])
  n <- nrow(x)
  stopifnot(ncol(x) == length(values))

  y <- as.matrix(x)
  xtx <- c(tcrossprod(values^2, y))
  xt1 <- c(tcrossprod(values, y))
  xt12 <- xt1^2

  cov_ <- \(x) {
    n <- nrow(x)
    stats::cov(x) * (n - 1) / n
  }

  theta <- cov_(cbind(xt1, xt12, xtx))

  grad_fun <- \(x) {
    a <- x[1]
    b <- x[2]
    c <- x[3]
    const <- (c - a^2 / r)^(-1)
    const * (r - 1)^(-1) * c(
      2 * a * ((b - a^2) * const / r - 1),
      1,
      -const * (b - a^2)
    )
  }

  grad <- grad_fun(c(mean(xt1), mean(xt12), mean(xtx)))
  c(t(grad) %*% theta %*% grad)
}
