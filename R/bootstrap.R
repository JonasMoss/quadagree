#' @keywords internal
ci_boot <- function(calc,
                    fun,
                    type,
                    transformer,
                    quants,
                    n_reps) {
  est_var <- fun(calc)
  est_t <- transformer$est(est_var$est)
  sd_t <- transformer$sd(est_var$est, sqrt(est_var$var))
  boots <- studentized_boots(
    calc,
    est_t,
    fun,
    transformer,
    n_reps
  )
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
studentized_boots <- function(calc,
                              est_t,
                              fun,
                              transformer,
                              n_reps) {
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
