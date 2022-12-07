library("future.apply")
plan(multisession)

## Parameters fro the simulation.
n <- 100
n_reps <- 10**4
r <- 5
mu <- rep(1, 5)
sigma <- matrix(0.5, r, r)
sigma <- cor(mtcars[, seq(r)])
diag(sigma) <- 1
p <- rep(1, 5)
p <- seq(1, 0.5, length.out = r)
mu <- 1:5
par <- conger_pop(mu, sigma)
f <- \(x) {
  ci = congerci::congerci(x)
  (c(ci[1] <= par) & (par <= ci[2]))
}

## The simulation itself.
#set.seed(313)
results <- simulate(n, mu, sigma, p, n_reps, f = f, laplace = FALSE)
mean(results)

x <- simulate(10000, mu, sigma, p, n_reps = 1, f = \(x) x, laplace = FALSE)[, , 1]
fleiss(x)
mat <- cov_mat(p, mu, sigma, gamma_est(NULL, sigma, type = "normal"))
avar_(mat, mu, sigma, TRUE)

x <- simulate(10000, mu, sigma, p, n_reps = 1, f = \(x) x, laplace = FALSE)[, , 1]
congerci::fleissci(x, type = "normal")


##
sqrt(avar_(mat, mu, sigma, TRUE))


## Parameters for the simulation.
n <- 10**4
n_reps <- 10**4
r <- 5
mu <- rep(1, 5)
mu <- 1:5
sigma <- matrix(0.3, r, r)
diag(sigma) <- 1
p <- rep(1, 5)

f <- function(x) {
  covariance <- cov(x, use = "pairwise.complete.obs")
  tr <- sum(diag(covariance))
  means <- colMeans(x, na.rm = TRUE)
  theta <- mean(means^2, na.rm = TRUE) - mean(means, na.rm = TRUE)^2
  c(sum(covariance) - tr, tr, theta)
}

## The simulation itself.
#set.seed(313)
results <- simulate(n, mu, sigma, p, n_reps, f = f, laplace = FALSE)
cov(t(results)) * n
cov_mat(p, mu, sigma, gamma_est(NULL, sigma, type = "normal"))





## The simulation itself.
#set.seed(313)
results <- simulate(n, mu, sigma, p, n_reps, f = fleiss, laplace = FALSE)
var(results) * n
mat <- cov_mat(p, mu, sigma, gamma_est(NULL, sigma, type = "normal"))
avar_(mat, mu, sigma, TRUE)

## The simulation itself.
#set.seed(313)
f <- \(x) c(congerci::fleiss(x), congerci:::avar(x, type = "normal", fleiss = TRUE))
results <- simulate(n, mu, sigma, p, n_reps, f = f, laplace = FALSE)
cov(t(results)) * n

mat <- cov_mat(p, mu, sigma, gamma_est(NULL, sigma, type = "normal"))
avar_(mat, mu, sigma, TRUE)









x <- simulate(100, mu, sigma, p, n_reps = 1, f = \(x) x, laplace = FALSE)[, , 1]
avar(x, type = "normal", fleiss = TRUE)


## Exact ci

## Parameters for the simulation.
n <- 10**4
n_reps <- 10**5
r <- 5
mu <- rep(1, 5)
#mu <- 1:5
sigma <- matrix(0.7, r, r)
diag(sigma) <- 1
p <- rep(1, 5)
std <- sqrt(avar_(mat, mu, sigma, TRUE))
par <- fleiss_pop(mu, sigma)
f <- function(x) {
  pm = 1.96 / sqrt(nrow(x)) * std
  (par >= fleiss(x) - pm) & (par <= fleiss(x) + pm)
}

## The simulation itself.
#set.seed(313)
results <- simulate(n, mu, sigma, p, n_reps, f = f, laplace = FALSE)
mean(results)
