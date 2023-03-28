x <- fleissci::dat.zapf2016[seq(20), 1]
y <- fleissci::dat.zapf2016[seq(20), 2]
n <- length(x)


indices <- arrangements::combinations(20, 10)
covs <- rep(NA, nrow(indices))

for(i in seq(nrow(indices))) {
  x_indices <- indices[i, ]
  y_indices <- setdiff(seq(n), indices[i, ])
  zx <- c(x[x_indices], y[y_indices])
  zy <- c(y[x_indices], x[y_indices])
  covs[i] = cov(zx, zy)
}


cor(x, y)
cor(c(x, y), c(y, x))
