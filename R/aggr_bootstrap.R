#' Studentized bootstrap estimates using transformers.
#'
#' @param n_reps Number of bootstrap repetitions.
#' @param x Data to estimate kappa on.
#' @param values Values for the different categories.
#' @param transformer A `transformer` object.
#' @param fleiss If `TRUE`, calculates Fleiss' kappa. If not, calculates
#'    Conger's kappa.
#' @return Studentized bootstrap estimates.
#' @keywords internal
studentized_boots_aggr <- function(n_reps,
                                   est_t,
                                   calc,
                                   transformer,
                                   est_fun,
                                   var_fun) {
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

#' @keywords internal
#' @rdname ci
ci_boot_aggr <- function(est,
                         sd,
                         calc,
                         transformer,
                         quants,
                         n_reps,
                         est_fun,
                         var_fun) {
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sd)
  boots <- studentized_boots_aggr(
    n_reps,
    est_t,
    calc,
    transformer,
    est_fun,
    var_fun
  )
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}
