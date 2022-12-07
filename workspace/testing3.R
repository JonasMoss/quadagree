library("future.apply")
plan(multisession)


## Parameters fro the simulation.
n <- 10**3
n_reps <- 10**5
r <- 5
mu <- seq(r)
sigma <- matrix(0.3, r, r)
diag(sigma) <- 1
p <- seq(1, 0.5, length.out = r)
f <- function(x) c(congerci::fleiss(x), congerci::conger(x))

## The simulation itself.
set.seed(313)
results <- simulate(n, mu, sigma, p, n_reps, f = f, laplace = FALSE)
cov(t(results)) * n

x <- simulate(10, mu, sigma, p, n_reps = 2, f = f, laplace = FALSE)
avar(x, type = "adf", fleiss = TRUE)
avar(x, type = "adf", fleiss = FALSE)

f <- function(x) {}
## The simulation itself.
set.seed(313)
results <- simulate(n, mu, sigma, p, n_reps, f = f, laplace = FALSE)
