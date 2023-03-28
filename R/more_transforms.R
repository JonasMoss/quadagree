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

  f <- \(i, j) {
    indices = c(y1 == i & y2 == j) | c(y1 == j & y2 == i) | c(y1 == i & y2 == i) | c(y1 == j & y2 == j)
    c(cor(y1[indices] == i, y2[indices] == i), var(y1[indices] == 1))
  }

  table(y1,)



  # w0 <- sapply(seq_along(y1), \(i) (y1[i] != y2[i]))
  # w0_org <- sapply(seq_along(y1), \(i) abs(y1[i] - y2[i])^2)
  # w <- w0/w0_org
  # non_na <- rep(NA, length(y1))
  # non_na[is.na(w)] = 0
  # non_na[!is.na(w)] = w[!is.na(w)]

  list(
    y1 = y1,
    y2 = y2,
    r1 = index1,
    r2 = index2,
    w = non_na
  )
}



data <- fleissci::dat.zapf2016
yy <- construct_y(data)
