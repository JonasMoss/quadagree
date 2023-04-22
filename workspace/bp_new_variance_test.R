### ============================================================================
###  Testing the second variance formulation in the document.
### ============================================================================

x = quadagree::dat.fleiss1971
values <- seq(ncol(x))

r <- sum(x[1, ])
y <- as.matrix(x)
calc <- bp_calc(y, values)
c1 <- bp_get_c1(y, values, type = 1)

bp_var_matrix1 <- \(calc, c1, r) {
  n <- nrow(calc)
  theta <- stats::cov(calc) * (n - 1) / n
  grad <- 1 / c1 * 2 / (r - 1) * c(-1, 1 / r)
  c(t(grad) %*% theta %*% grad)
}

bp_var_matrix2 <- \(calc, c1, r) {
  n <- nrow(calc)
  phi <- stats::cov(calc) * (n - 1) / n
  # Modify the elements of phi to make the variance calculation easier.
  phi[1, 1] = phi[1, 1]
  phi[1, 2] = phi[2, 1] = -phi[1, 2] / r
  phi[2, 2] = phi[2, 2] / r^2
  1 / c1 ^ 2 * 4 / (r - 1) ^ 2 * sum(phi)
}

bp_var_matrix1(calc, c1, r)
bp_var_matrix2(calc, c1, r)

microbenchmark::microbenchmark(bp_var_matrix1(calc, c1, r),
                               bp_var_matrix2(calc, c1, r))
