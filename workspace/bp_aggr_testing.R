### ============================================================================
###   Testing inference for aggregated Brennan-Prediger.
### ============================================================================

# Testing data
x <- fleissci::dat.fleiss1971

#' Convert data on raw form to Fleiss form.
#'
#' @param x Data on raw form.
#' @param n_cat Number of categories.
#' @return Data on Fleiss form.
raw_to_fleiss <- \(x, n_cat) {
  n <- nrow(x)
  y <- matrix(0, nrow = n, ncol = n_cat)
  for (i in seq(n)) {
    for (j in x[i, ]) {
      y[i, j] = 1 + y[i, j]
    }
  }
  y
}


bp1 <- \(x, values = seq(ncol(x))) {
  r <- sum(x[1, ])
  n_cat <- ncol(x)
  w <- outer(values, values, Vectorize(\(x, y) (x - y)^2))
  c1 = sum(w) / n_cat ^ 2

  y <- as.matrix(x)
  xtx <- c(tcrossprod(values ^ 2, y))
  xt12 <- c(tcrossprod(values, y)) ^ 2

  sxtx <- mean(xtx)
  sxt12 <- mean(xt12)

  disagreement <- 2/(r - 1) * (sxtx - 1/r * sxt12)
  1 - disagreement / c1
}

bp1_var <- \(x, values = seq(ncol(x))) {
  r <- sum(x[1, ])
  n_cat <- ncol(x)
  n <- nrow(x)
  w <- outer(values, values, Vectorize(\(x, y) (x - y)^2))
  c1 = sum(w) / n_cat ^ 2

  y <- as.matrix(x)
  xtx <- c(tcrossprod(values ^ 2, y))
  xt12 <- c(tcrossprod(values, y)) ^ 2

  theta <- cov(cbind(xtx, xt12)) * (n - 1) / n

  grad <- 1 / c1 * 2 / (r - 1) * c(-1, 1/r)
  c(t(grad) %*% theta %*% grad)
}

true_dist <- c(0.1,0.1,0.5,0.3)
skills <- c(0.9, 0.9, 0.8)

x <- agreeable::simulate_jsm(50000, skills, model = "fleiss",
                             true_dist = true_dist)
x <- raw_to_fleiss(x, n_cat = length(true_dist))
bp1_var(x, 1:4)

results <- replicate(20000, {
  x <- agreeable::simulate_jsm(200, skills, model = "fleiss",
                               true_dist = true_dist)
  x <- raw_to_fleiss(x, n_cat = length(true_dist))
  bp1(x, 1:4)
})

var(results) * 200
