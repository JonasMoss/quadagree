x <- dat.fleiss1971
x <- rbind(x, x, x, x, x, x, x, x, x, x)
x <- rbind(x, x, x, x, x, x, x, x, x, x)
microbenchmark::microbenchmark(fleissci_aggr(x))

microbenchmark::microbenchmark({
  colSums(grad * (theta %*% grad))
}, {
  c(crossprod(grad, theta %*% grad))
}, {
  c(crossprod(grad, crossprod(theta, grad)))
}, {
  c(tcrossprod(theta %*% grad, grad))
}, {
  t(grad) %*% theta %*% grad
})

n <- nrow(x)
microbenchmark::microbenchmark(
  dplyr::sample_n(x, n, replace = TRUE),
  x[sample(n, replace = TRUE), ]
)
