bp_prepare <- \(x, values, type) {
  r <- sum(x[1, ])
  y <- as.matrix(x)
  calc <- bp_calc(y, values)
  c1 <- bp_get_c1(y, values, type)
  list(calc = calc, c1 = c1, r = r)
}

bp_calc <- \(y, values) {
  xtx <- c(tcrossprod(values^2, y))
  xt12 <- c(tcrossprod(values, y))^2
  cbind(xtx, xt12)
}

bp_est <- \(x, values, type) {
  args <- bp_prepare(x, values, type)
  do.call(bp_est_matrix, args)
}

bp_est_matrix <- \(calc, c1, r) {
  means = colMeans(calc)
  disagreement <- 2 / (r - 1) * (means[1] - 1 / r * means[2])
  unname(1 - disagreement / c1)
}

bp_get_c1 <- \(x, values, type) {
  if (type == 1) {
    w <- outer(values, values, Vectorize(\(x, y) (x - y)^2))
    n_cat <- ncol(x)
    sum(w) / n_cat^2
  } else {
    0.5 * (max(values) - min(values))^2
  }
}

bp_var <- \(x, values, type) {
  args <- bp_prepare(x, values, type)
  do.call(bp_var_matrix, args)
}

bp_var_matrix <- \(calc, c1, r) {
  n <- nrow(calc)
  theta <- stats::cov(calc) * (n - 1) / n
  grad <- 1 / c1 * 2 / (r - 1) * c(-1, 1 / r)
  c(t(grad) %*% theta %*% grad)
}
