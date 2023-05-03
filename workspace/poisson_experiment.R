res <- replicate(1000, length(unique(sample.int(n = 100, replace = TRUE))))

hist(res)

res <- replicate(1000, sum(rpois(100, lambda = 1)))
mean(res)

n <- 100
microbenchmark::microbenchmark(
  sample.int(n, replace = TRUE),
  rpois(n, 1)
  )



set.seed(313)
n <- 1000
x <- rnorm(n)




microbenchmark::microbenchmark({
  res1 <- replicate(100, {
    indices <- rpois(n, lambda = 1)
    sum(x * indices) / length(indices)
  })
}, {
  res2 <- replicate(100, {
    indices <- sample.int(n = n, replace = TRUE)
    sum(x[indices]) / n
  })
})


quantile(res1, c(0.025, 0.975), na.rm = TRUE)
quantile(res2, c(0.025, 0.975), na.rm = TRUE)
t.test(x)$conf.int


z <- psych::bfi[, 1:5]
colSums(t(z) * c(1, 10, 100, 1000, 10000))
Rfast::Table(rowSums(t(t(z) * c(1, 10, 100, 1000, 10000))))
