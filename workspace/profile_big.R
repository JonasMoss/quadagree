x <- dat.zapf2016
x <- rbind(x, x, x, x, x, x, x, x, x, x)
x <- rbind(x, x, x, x, x, x, x, x, x, x)
nrow(x) #5000

profvis::profvis(fleissci(x, boots = TRUE))

x <- dat.fleiss1971
x <- rbind(x, x, x, x, x, x, x, x, x, x)
x <- rbind(x, x, x, x, x, x, x, x, x, x)
x <- rbind(x, x, x, x, x, x, x, x, x, x)
profvis::profvis(fleissci_aggr(x, boot = TRUE))

microbenchmark::microbenchmark(fleissci_aggr(x), times = 10000)

#microbenchmark::microbenchmark(fleissci(x), times = 10000)

x <- dat.zapf2016
i_row <- \(n) unlist(lapply(seq_len(n), seq.int, n))
i_col <- \(n) rep.int(seq_len(n), times = rev(seq_len(n)))
rows <- i_row(ncol(x))
cols <- i_col(ncol(x))
y <- t(x) - colMeans(x, na.rm = TRUE)
z <- y[cols, ] * y[rows, ]


### ============================================================================
###  Note.
###
###  It could be possible to speed up the bootstrap by precomputing
###    t(x)[cols, ] * t(x)[rows, ].
###  Maybe we have to use tricks to avoid excessive copying. This can be done
###   using environments.
### ============================================================================

y <- t(x)

x_cols <- x[, cols]
x_rows <- x[, rows]
z <- x_cols * x_rows

x_cols <- y[cols, ]
x_rows <- y[rows, ]
z <- x_cols * x_rows



microbenchmark::microbenchmark(
  {
    mat <- z - rowMeans(z, na.rm = TRUE)
    tcrossprod(mat)
  },
  cov(y)
)

x <- rnorm(10000)
microbenchmark::microbenchmark(
  list(mean(x), var(x)),
  {
  n <- length(x)
  mu <- sum(x) / n
  list(mu, sum((x - mu)^2) / n)
  },
  var(x),
  sum(x) / n,
  times = 1000
)

bench::bench_time( {
  n <- length(x)
  mu <- sum(x) / n
  list(mu, sum((x - mu)^2) / (n - 1))
})

x <- rlnorm(10000)
bench::bench_time( {
  mean(log(x))
})

bench::bench_time( {
  mean(x)
})


r = 10
r ^ 2 * (mean(mu^2) - mean(mu)^2)

sum(apply(arrangements::combinations(mu, 2), 1, \(mu) (mu[2] - mu[1])^2))



zz <- fleiss_aggr_prepare(dat.fleiss1971, 1:5)
zz <- zz$xx
mu <- colMeans(zz)
microbenchmark::microbenchmark(
  cov(zz),
  tcrossprod(t(zz) - mu) / nrow(zz),
  times = 1000
)

cov(zz)
tcrossprod(t(zz) - mu) / nrow(zz)
