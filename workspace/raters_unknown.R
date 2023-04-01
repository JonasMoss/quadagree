tr <- \(x) {
  sum(diag(x))
}

cov_ <- \(x) {
  n <- nrow(x)
  cov(x) * (n - 1) / n
}

var_ <- \(x) {
  n <- length(x)
  var(x) * (n - 1) / n
}

#' Quadratic Fleiss' kappa for data on Fleiss form.
#'
#' @param x Data on Fleiss form.
#' @param values Values to attach to each column on the Fleiss form data.
#'    Defaults to `1:C`, where `C` is the number of categories.
#' @return Calculated value of Fleiss' kappa.
fleiss_aggr <- \(x, values = seq(ncol(x))) {
  r <- sum(x[1, ])
  stopifnot(ncol(x) == length(values))

  xtx <- apply(x, 1, \(row) sum(values ^ 2 * row))
  xt1 <- apply(x, 1, \(row) sum(values * row))

  vxt1 <- var_(xt1)
  extx <- mean(xtx)
  ext1 <- mean(xt1)
  ext2 <- mean(xt1^2)

  1 / (r - 1) * (vxt1 / (extx - ext1^2 / r) - 1)
  1 / (r - 1) * ((ext2 - ext1^2) / (extx - ext1^2 / r) - 1)
}

fleiss_aggrci <- \(x, values = seq(ncol(x))) {
  r <- sum(x[1, ])
  n <- nrow(x)
  stopifnot(ncol(x) == length(values))

  xtx <- apply(x, 1, \(row) sum(values ^ 2 * row))
  xt1 <- apply(x, 1, \(row) sum(values * row))
  xt12 <- xt1^2

  theta <- cov_(cbind(xt1, xt12, xtx))

  grad_fun <- \(x) {
    a <- x[1]
    b <- x[2]
    c <- x[3]
    const <- (c - a^2/r)^(-1)

    const * (r-1)^(-1)*c(
      2 * a * ((b - a^2) * const / r - 1),
      1,
      -const * (b - a^2)
    )
  }

  grad <- grad_fun(c(mean(xt1), mean(xt12), mean(xtx)))
  c(t(grad) %*% theta %*% grad)
}


## Testing

long_to_wide <- function(x, rating = "rating") {
  as.matrix(dplyr::select(tidyr::spread(x, key = "judge", value = rating), -1))
}

long_to_wide(fleissci::dat.zapf2016)















r <- 6
fun <- \(x) {
  a <- x[1]
  b <- x[2]
  c <- x[3]
  1 / (r - 1) * ((b - a^2) / (c - a^2 / r) - 1)
}

grad_fun <- \(x) {
  a <- x[1]
  b <- x[2]
  c <- x[3]
  const <- (c - a^2/r)^(-1)

  const * (r-1)^(-1)*c(
    2 * a * ((b - a^2) * const / r - 1),
    1,
    -const * (b - a^2)
  )
}
numDeriv::grad(fun, c(0.5, 0.6, 0.7))
grad_fun(c(0.5, 0.6, 0.7))



x <- fleissci::dat.fleiss1971

#x <- t(t(xx)*(1:6))

xtx <- apply(x, 1, \(row) sum((1:5)^2 * row))
xt1 <- apply(x, 1, \(row) sum((1:5) * row))



xt1 <- rowSums(x)

var1tx <- var(xt1)

c_xtx_xt1 <- mean(xtx * xt1) - (tr(sigma) + mutmu)









z <- fake_wide(fleissci::dat.fleiss1971)
z_ <- z
for(i in 1:nrow(z)) {
  z_[i, 1:6] <- z[i, sample(1:6, 6)]
}

tr(cov_(z)) + crossprod(colMeans(z))
tr(cov_(z_)) + crossprod(colMeans(z_))
sum(cov_(z))
sum(cov_(z_))

mean(colMeans(z_))
mean(colMeans(z))
