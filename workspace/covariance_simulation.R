### ===========================================================================
###   Covariance of residuals simulation using dat.zapf2016
### ===========================================================================

data <- fleissci::dat.zapf2016

sampler <- \(n) dplyr::slice_sample(data, n = n)

construct_y <- \(data, vars = NULL) {
  r <- ncol(data)
  n <- nrow(data)
  indices <- arrangements::combinations(r, 2)
  indices <- rbind(indices, cbind(indices[, 2], indices[, 1]))
  y1 <- y2 <- c()
  index1 <- index2 <- c()
  for (row in seq(nrow(indices))) {
    y1 <- c(y1, data[, indices[row, 1]])
    y2 <- c(y2, data[, indices[row, 2]])
    index1 <- c(index1, rep(indices[row, 1], n))
    index2 <- c(index2, rep(indices[row, 2], n))
  }
  list(
    y1 = y1,
    y2 = y2,
    r1 = index1,
    r2 = index2,
    i = rep(seq(n), r * (r - 1))
  )
}

yy <- construct_y(data)
n <- nrow(data)
r <- ncol(data)
residuals <- lm(yy$y1 ~ yy$y2)$residuals
resids = list()
for (x in unique(yy$y1)) {
  resids[[x]] = residuals[yy$y2 == x]
}

