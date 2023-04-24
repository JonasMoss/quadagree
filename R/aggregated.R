bp_aggr_prepare <- \(x, values, kind) {
  r <- sum(x[1, ])
  y <- as.matrix(x)
  c1 <- bp_aggr_get_c1(values, kind)
  xtx <- c(tcrossprod(values^2, y))
  xt12 <- c(tcrossprod(values, y))^2
  list(
    xx = cbind(xtx, xt12),
    n = nrow(x),
    r = r,
    c1 = c1
  )
}

bp_aggr_get_c1 <- \(values, kind) {
  if (kind == 1) {
    w <- outer(values, values, Vectorize(\(x, y) (x - y)^2))
    n_cat <- length(values)
    sum(w) / n_cat^2
  } else {
    0.5 * (max(values) - min(values))^2
  }
}

bp_aggr_est_matrix <- \(calc) {
  means <- colMeans(calc$xx)
  disagreement <- 2 / (calc$r - 1) * (means[1] - 1 / calc$r * means[2])
  unname(1 - disagreement / calc$c1)
}

bp_aggr_var_matrix <- \(calc) {
  phi <- stats::cov(calc$xx) * (calc$n - 1) / calc$n
  phi[1, 1] <- phi[1, 1]
  phi[1, 2] <- phi[2, 1] <- -phi[1, 2] / calc$r
  phi[2, 2] <- phi[2, 2] / calc$r^2
  1 / calc$c1^2 * 4 / (calc$r - 1)^2 * sum(phi)
}

fleiss_aggr_prepare <- \(x, values) {
  r <- sum(x[1, ])
  y <- as.matrix(x)
  xtx <- c(tcrossprod(values^2, y))
  xt1 <- c(tcrossprod(values, y))
  xt12 <- xt1^2

  list(
    xx = cbind(xt1, xt12, xtx),
    n = nrow(x),
    r = r
  )
}

fleiss_aggr_var <- \(calc) {
  theta <- stats::cov(calc$xx) * (calc$n - 1) / calc$n
  a <- mean(calc$xx[, 1])
  b <- mean(calc$xx[, 2])
  c <- mean(calc$xx[, 3])
  k <- (c - a^2 / calc$r)^(-1)
  grad <- k / (calc$r - 1) *
    c(2 * a * ((b - a^2) * k / calc$r - 1), 1, -k * (b - a^2))
  c(crossprod(grad, theta %*% grad))
}

fleiss_aggr_est <- \(calc) {
  ext1 <- mean(calc$xx[, 1])
  ext2 <- mean(calc$xx[, 2])
  extx <- mean(calc$xx[, 3])
  1 / (calc$r - 1) * ((ext2 - ext1^2) / (extx - ext1^2 / calc$r) - 1)
}
