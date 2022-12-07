## Parameters fro the simulation.
n <- 10**5
n_reps <- 10**5
r <- 3
mu <- seq(r)
sigma <- matrix(c(
  1, 0.5, 0.6,
  0.5, 1, 0.7,
  0.6, 0.7, 1
), nrow = 3)

p <- seq(1, 0.5, length.out = r)

x <- simulate(n, mu, sigma, p, 1, f = \(x) x)[, , 1]
n <- nrow(x)


mu <- colMeans(x, na.rm = TRUE)
sigma <- cov(x, use = "pairwise.complete.obs") * (n - 1) / n

gamma <- gamma_est(x, sigma)
gamma_est(x, sigma, type = "normal")
gamma_est(x, sigma, type = "elliptical")
