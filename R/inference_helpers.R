#' @keywords internal
quadagree_internal <- function(calc,
                               transform,
                               conf_level,
                               alternative = c("two.sided", "greater", "less"),
                               bootstrap,
                               n_reps,
                               call,
                               fun,
                               ...) {
  alternative <- match.arg(alternative)
  transformer <- get_transformer(transform)
  quants <- limits(alternative, conf_level)

  calc <- list(xx = x, n = nrow(x))
  est_var <- fun(calc)
  est <- est_var$est
  sd <- sqrt(est_var$var)

  ci <- if (!bootstrap) {
    ci_asymptotic(est, sd, nrow(x), transformer, quants)
  } else {
    ci_boot(calc, fun, transformer, quants, n_reps)
  }

  names(ci) <- quants
  attr(ci, "conf_level") <- conf_level
  attr(ci, "alternative") <- alternative
  attr(ci, "type") <- type
  attr(ci, "n") <- nrow(x)
  attr(ci, "transform") <- transform
  attr(ci, "bootstrap") <- bootstrap
  attr(ci, "n_reps") <- n_reps
  attr(ci, "estimate") <- est
  attr(ci, "sd") <- sd
  attr(ci, "call") <- call
  class(ci) <- "quadagree"
  ci[2] <- min(ci[2], 1)
  ci
}

#' @keywords internal
ci_asymptotic <- function(est, sd, n, transformer, quants) {
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sd)
  multiplier <- stats::qt(quants, n - 1) / sqrt(n - 1)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
ci_boot <- function(calc, fun, transformer, quants, n_reps) {
  est_var <- fun(calc)
  est_t <- transformer$est(est_var$est)
  sd_t <- transformer$sd(est_var$est, sqrt(est_var$var))
  boots <- bootstrapper(calc, est_t, fun, transformer, n_reps)
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
bootstrapper <- function(calc, est_t, fun, transformer, n_reps) {
  calc_new <- calc
  trans_est <- transformer$est
  trans_sd <- transformer$sd
  future.apply::future_replicate(n_reps,
    {
      indices_star <- sample.int(calc$n, replace = TRUE)
      calc_new$xx <- calc$xx[indices_star, ]
      est_var <- fun(calc_new)
      est_star <- est_var$est
      sd_star <- sqrt(est_var$var)
      (trans_est(est_star) - est_t) / trans_sd(est_star, sd_star)
    },
    future.seed = TRUE
  )
}
