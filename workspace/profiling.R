#' Calculates the empirical capital pi matrix.
#'
#' @param x Vector of probabilities for being missing.
#' @return The capital pi matrix.
#' @keywords internal
pi_mat_empirical <- \(x) {
  r <- ncol(x)
  ind2 <- arrangements::combinations(r, 2, replace = TRUE)
  ind4 <- arrangements::combinations(seq_len(nrow(ind2)), 2, replace = TRUE)
  p2_hats <- apply(ind2, 1, \(i) mean(!is.na(x[, i[1]]) & !is.na(x[, i[2]])))
  hats <- apply(ind4, 1, \(i) {
    j1 <- ind2[i[1], ]
    j2 <- ind2[i[2], ]
    mean(!is.na(x[, j1[1]]) & !is.na(x[, j1[2]]) &
           !is.na(x[, j2[1]]) & !is.na(x[, j2[2]])) /
      (p2_hats[i[1]] * p2_hats[i[2]])
  })

  new_mat <- matrix(NA, choose(r + 1, 2), choose(r + 1, 2))
  new_mat[ind4] <- hats
  as.matrix(Matrix::forceSymmetric(new_mat, uplo = "U"))
}

#' Calculates the empirical capital pi matrix.
#'
#' @param x Vector of probabilities for being missing.
#' @return The capital pi matrix.
#' @keywords internal
pi_mat_empirical2 <- \(x) {
  r <- ncol(x)
  ind2 <- arrangements::combinations(r, 2, replace = TRUE)
  ind4 <- arrangements::combinations(seq_len(nrow(ind2)), 2, replace = TRUE)

  nisna <- !is.na(x)

  p2_hats <- apply(ind2, 1, \(i) mean(nisna[, i[1]] & nisna[, i[2]]))

  hats <- apply(ind4, 1, \(i) {
    j1 <- ind2[i[1], ]
    j2 <- ind2[i[2], ]
    mean(nisna[, j1[1]] & nisna[, j1[2]] &
           nisna[, j2[1]] & nisna[, j2[2]]) /
      (p2_hats[i[1]] * p2_hats[i[2]])
  })

  new_mat <- matrix(NA, choose(r + 1, 2), choose(r + 1, 2))
  new_mat[ind4] <- hats
  as.matrix(Matrix::forceSymmetric(new_mat, uplo = "U"))
}

pi_mat_empirical3 <- \(x) {
  r <- ncol(x)
  if(!anyNA(x)) {
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



x <- dat.zapf2016
all.equal(pi_mat_empirical(x), pi_mat_empirical2(x))
all.equal(pi_mat_empirical(x), pi_mat_empirical3(x))
microbenchmark::microbenchmark(
  pi_mat_empirical(x),
  pi_mat_empirical2(x),
  pi_mat_empirical3(x)
)



get_diag_indices <- function(r) {
  e_mat <- matrixcalc::elimination.matrix(r)
  indices <- rep(0, r^2)
  indices[c(1, (r + 1) * (1:(r - 1)) + 1)] <- 1
}

get_diag_indices2 <- function(r) {
  indices <- rep(0, choose(r + 1, 2))
  indices[c(1, 1 + cumsum(r:2))] <- 1
  indices
}

microbenchmark::microbenchmark(
  get_diag_indices(5),
  get_diag_indices2(5),
  get_diag_indices(10),
  get_diag_indices2(10)
)

all.equal(get_diag_indices(5), get_diag_indices2(5))
all.equal(get_diag_indices(10), get_diag_indices2(10))

### Trace


tr <- \(x) {
  sum(x[1L + 0L:(dim(x)[1L] - 1L) * (dim(x)[1L] + 1L)])
}

tr <- \(x) sum(diag(x))
tr2 <- \(x) sum(x) - 2 * sum(Rfast::lower_tri(x))
tr3 <- \(x) sum(x) - 2 * sum(Rfast::lower_tri(x, suma = TRUE))
x <-  cov(dat.zapf2016)

tr(x)
tr2(x)
tr3(x)
tr4(x)
microbenchmark::microbenchmark(
  tr(x),
  tr2(x),
  tr3(x),
  tr4(x),
  times = 100
)

y <- as.matrix(x)
microbenchmark::microbenchmark(
  colMeans(y),
  Rfast::colmeans(y)
)

mu <- colMeans(y)
n <- length(mu)
microbenchmark::microbenchmark(
  sum(mu^2) / n,
  sum(mu)^2 / n^2,
  sum(diag(y)),
  sum(y)
)
mean(mu^2) - mean(mu)^2

x <- as.matrix(dat.zapf2016)
microbenchmark::microbenchmark(
  colSums(x) / nrow(x),
  colMeans(x)
)
