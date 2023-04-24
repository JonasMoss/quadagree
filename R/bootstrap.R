#' @keywords internal
ci_boot <- function(calc,
                    est_fun,
                    var_fun,
                    type,
                    transformer,
                    quants,
                    n_reps) {
  est <- est_fun(calc)
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sqrt(var_fun(calc)))
  boots <- studentized_boots(
    calc,
    est_t,
    est_fun,
    var_fun,
    transformer,
    n_reps
  )
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
studentized_boots <- function(calc,
                              est_t,
                              est_fun,
                              var_fun,
                              transformer,
                              n_reps) {
  calc_new <- calc
  trans_est <- transformer$est
  trans_sd <- transformer$sd
  future.apply::future_replicate(n_reps,
    {
      indices_star <- sample(calc$n, replace = TRUE)
      calc_new$xx <- calc$xx[indices_star, ]
      est_star <- est_fun(calc_new)
      sd_star <- sqrt(var_fun(calc_new))
      (trans_est(est_star) - est_t) / trans_sd(est_star, sd_star)
    },
    future.seed = TRUE
  )
}
