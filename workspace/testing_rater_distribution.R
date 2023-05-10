sigma <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.7, 0.3, 0.7, 1), 3, 3)
mu <- c(1, 2, 3)

bp_pop(mu, sigma, c(1, 3), 2)



x <- dat.zapf2016
n <- length(x)
r <- ncol(x)
means <- colMeans(x)
ind <- arrangements::combinations(r, 2)

sum(apply(ind, 1, \(x) (means[x[1]] - means[x[2]])^2)) / r^2
mean(means^2) - mean(means)^2
