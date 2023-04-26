bp_aggr_prepare <- \(x, values, kind) {
  y <- as.matrix(x)
  r <- sum(y[1, ])
  c1 <- bp_get_c1(values, kind)
  xtx <- c(tcrossprod(values^2, y))
  xt12 <- c(tcrossprod(values, y))^2
  list(
    xx = cbind(xtx, xt12),
    n = nrow(x),
    r = r,
    c1 = c1
  )
}

bp_aggr_fun <- \(calc) {
  means <- colMeans(calc$xx)
  calcr_inv_1 <- 1 / (calc$r - 1)
  calcr_inv <- 1 / calc$r
  calcc1_inv <- 1 / calc$c1
  disagreement <- 2 * calcr_inv_1 * (means[1] - calcr_inv * means[2])
  phi <- stats::cov(calc$xx) * (calc$n - 1) / calc$n
  k <- phi[1, 1] - 2 * phi[1, 2] * calcr_inv + phi[2, 2] * calcr_inv^2

  est <- unname(1 - disagreement * calcc1_inv)
  var <- calcc1_inv^2 * 4 * calcr_inv_1^2 * k
  list(est = est, var = var)
}

fleiss_aggr_prepare <- \(x, values) {
  y <- as.matrix(x)
  r <- sum(y[1, ])
  xtx <- c(tcrossprod(values^2, y))
  xt1 <- c(tcrossprod(values, y))
  xt12 <- xt1^2

  list(
    xx = cbind(xt1, xt12, xtx),
    n = nrow(x),
    r = r
  )
}

fleiss_aggr_fun <- \(calc) {
  theta <- stats::cov(calc$xx) * (calc$n - 1) / calc$n
  m <- colMeans(calc$xx)

  k <- calc$r / (m[3] * calc$r - m[1]^2)
  km <- k * (m[2] - m[1]^2)
  calc_r_inv <- 1 / (calc$r - 1)

  grad <- k * calc_r_inv * c(2 * m[1] * (km / calc$r - 1), 1, -km)

  est <- calc_r_inv * (km - 1)
  var <- c(crossprod(grad, theta %*% grad))

  list(est = est, var = var)
}
