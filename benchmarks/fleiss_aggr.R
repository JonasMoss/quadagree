library("fleissci")

x <- dat.fleiss1971
microbenchmark::microbenchmark(
  fleiss_aggr_var(x),
  fleiss_aggr_var2(x)
)

fleiss_aggr2 <- \(x, values = seq(ncol(x))) {
  r <- sum(x[1, ])
  stopifnot(ncol(x) == length(values))

  y <- as.matrix(x)
  xtx <- tcrossprod(values ^ 2, y)
  xt1 <- tcrossprod(values, y)

  extx <- mean(xtx)
  ext1 <- mean(xt1)
  ext2 <- mean(xt1^2)

  1 / (r - 1) * ((ext2 - ext1^2) / (extx - ext1^2 / r) - 1)
}

#' @return Calculated value of Fleiss' kappa.
fleiss_aggr_var2 <- \(x, values = seq(ncol(x))) {
  r <- sum(x[1, ])
  n <- nrow(x)
  stopifnot(ncol(x) == length(values))

  values2 <- values^2
  y <- as.matrix(x)
  xtx <- c(tcrossprod(values ^ 2, y))
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


values = seq(ncol(x))
r <- sum(x[1, ])
xtx <- apply(x, 1, \(row) sum(values^2 * row))
xt1 <- apply(x, 1, \(row) sum(values * row))

as.matrix(x) %*% values
