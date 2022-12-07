#' Simulate data from multivariate normal with missing observations.
#'
#' @param n Number of observations.
#' @param mu Mean vector.
#' @param sigma Covariance matrix.
#' @param p Vector of probabilities for being missing.
#' @param n_reps The number of repetitions desired.
#' @param f Function applied to each row.
#' @param laplace If `TRUE`, simulates from the multivariate Laplace.
#' @return Simulated values with `f` applied to them.
simulate <- function(n, mu, sigma, p, n_reps, f = mean, laplace = FALSE) {
  future.apply::future_replicate(n_reps,
    {
      z <- if (laplace) {
        LaplacesDemon::rmvl(n, mu = mu, Sigma = sigma)
      } else {
        MASS::mvrnorm(n, mu = mu, Sigma = sigma)
      }
      for (i in seq(r)) {
        indices <- sample(n, n - p[i] * n)
        z[indices, i] <- NA
      }
      f(z)
    },
    future.seed = TRUE
  )
}
