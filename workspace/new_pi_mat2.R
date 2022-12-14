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

pi_mat(c(1, 0.75, 0.5), vech = TRUE)


pi_mat2 <- \(p) {
  r <- length(p)
  f <- \(x) prod(p[unique(x)])
  indices <- arrangements::combinations(r, 2, replace = TRUE)
  ps <- apply(indices, 1, f)
  g <- Vectorize(\(x, y) f(c(indices[x, ], indices[y, ])) /  (ps[x] * ps[y]))
  x <- seq(nrow(indices))
  outer(x, x, g)
}

p = seq(1, 0.5, length.out = 7)
microbenchmark::microbenchmark(pi_mat(p, vech = TRUE), pi_mat2(p))
pi_mat2(c(1, 0.75, 0.5))

## Parameters fro the simulation.
n <- 10**4
n_reps <- 10**5
p = seq(1, 0.5, length.out = 3)
r = length(p)
mu <- seq(r)
sigma <- matrix(0.7, r, r)
x <- simulate(n, mu, sigma, p, n_reps = 1, f = \(x) x, laplace = FALSE)[, , 1]


pi_mat_empirical <- \(x) {
  r <- ncol(x)
  indices2 <- arrangements::combinations(r, 2, replace = TRUE)
  indices4 <- arrangements::combinations(seq(nrow(indices2)), 2, replace = TRUE)

  p2_hats <- apply(indices2, 1, \(i) mean(!is.na(x[, i[1]]) & !is.na(x[, i[2]])))
  hats <- apply(indices4, 1, \(i) {
    j1 = indices2[i[1], ]
    j2 = indices2[i[2], ]
    mean(!is.na(x[, j1[1]]) & !is.na(x[, j1[2]]) &
           !is.na(x[, j2[1]]) & !is.na(x[, j2[2]])) /
      (p2_hats[i[1]] * p2_hats[i[2]])
  })

  new_mat <- matrix(NA, choose(r+1, 2), choose(r+1, 2))
  new_mat[indices4] <- hats
  as.matrix(Matrix::forceSymmetric(new_mat,uplo="U"))

}



#p_hat = colMeans(!is.na(x))
p2_hat <- apply(indices2, 1, \(i) mean(!is.na(x[, i[1]]) & !is.na(x[, i[2]])))
p4_hat <- apply(indices4, 1, \(i) {
  j1 = indices2[i[1], ]
  j2 = indices2[i[2], ]
  mean(!is.na(x[, j1[1]]) & !is.na(x[, j1[2]]) &
         !is.na(x[, j2[1]]) & !is.na(x[, j2[2]]))
})

new_mat <- matrix(NA, choose(r+1, 2), choose(r+1, 2))
new_mat[indices4] <- p4_hat

