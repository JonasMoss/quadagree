#' Calculates asymptotic confidence intervals.
#'
#' @param x Data to estimate kappa on.
#' @param est,sd The estimate and estimated standard deviation.
#' @param n Number of observations.
#' @param type Type of confidence interval. Either `adf`, `elliptical`, or
#'   `normal`.
#' @param transformer A transformer object.
#' @param quants Quantiles for the confidence interval.
#' @param n_reps Number of bootstrap samples if `bootstrap = TRUE`. Ignored if
#'   `bootstrap = FALSE`.
#' @param fleiss If `TRUE`, calculates Fleiss' kappa. If not, calculates
#'    Conger's kappa.
#' @keywords internal
#' @name ci
ci_asymptotic <- function(est, sd, n, transformer, quants) {
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sd)
  multiplier <- stats::qt(quants, n - 1) / sqrt(n - 1)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' @keywords internal
#' @rdname ci
ci_boot <- function(x,
                    est,
                    sd,
                    type,
                    transformer,
                    quants,
                    n_reps,
                    fleiss = FALSE,
                    kind = NULL,
                    values = NULL) {
  boots <- studentized_boots(n_reps, x, type, transformer, fleiss, kind, values)
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sd)
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

#' Studentized bootstrap estimates using transformers.
#'
#' @param n_reps Number of bootstrap repetitions.
#' @param x Data to estimate kappa on.
#' @param type Type of confidence interval. Either `adf`, `elliptical`, or
#'   `normal`.
#' @param transformer A `transformer` object.
#' @param fleiss If `TRUE`, calculates Fleiss' kappa. If not, calculates
#'    Conger's kappa.
#' @return Studentized bootstrap estimates.
#' @keywords internal
studentized_boots <- function(n_reps,
                              x,
                              type,
                              transformer,
                              fleiss = FALSE,
                              kind = NULL,
                              values = NULL) {
  if (!is.null(kind)) {
    est_fun <- bp
    c1 <- bp_aggr_get_c1(values, kind)
    var_fun <- \(y) avar_bp(y, type, c1)
  } else {
    if (fleiss) {
      est_fun <- methods::getFunction("fleiss")
    } else {
      est_fun <- conger
    }
    var_fun <- \(y) avar(y, type, fleiss)
  }

  est <- est_fun(x)
  future.apply::future_replicate(n_reps,
    {
      indices_star <- sample(nrow(x), nrow(x), replace = TRUE)
      est_star <- est_fun(x[indices_star, ])
      sd_star <- sqrt(var_fun(x[indices_star, ]))
      (transformer$est(est_star) - transformer$est(est)) /
        transformer$sd(est_star, sd_star)
    },
    future.seed = TRUE
  )
}

#' @keywords internal
#' @rdname ci
ci_boot_aggr <- function(x,
                         values,
                         est,
                         sd,
                         transformer,
                         quants,
                         n_reps,
                         fleiss = FALSE) {
  boots <- studentized_boots_aggr(n_reps, x, values, transformer)
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sd)
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

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
                                   x,
                                   values,
                                   transformer) {
  est <- fleiss_aggr(x, values)
  future.apply::future_replicate(n_reps,
    {
      indices_star <- sample(nrow(x), nrow(x), replace = TRUE)
      est_star <- fleiss_aggr(x[indices_star, ], values)
      sd_star <- sqrt(fleiss_aggr_var(x[indices_star, ], values))
      (transformer$est(est_star) - transformer$est(est)) /
        transformer$sd(est_star, sd_star)
    },
    future.seed = TRUE
  )
}


#' @keywords internal
#' @rdname ci
ci_boot_bp_aggr <- function(x,
                            values,
                            kind,
                            est,
                            sd,
                            transformer,
                            quants,
                            n_reps) {
  est_t <- transformer$est(est)
  sd_t <- transformer$sd(est, sd)
  boots <- studentized_boots_bp_aggr(n_reps, est_t, x, values, kind, transformer)
  multiplier <- stats::quantile(boots, quants, na.rm = TRUE)
  sort(transformer$inv(est_t + multiplier * sd_t))
}

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
studentized_boots_bp_aggr <- function(n_reps,
                                      est_t,
                                      x,
                                      values,
                                      kind,
                                      transformer) {
  future.apply::future_replicate(n_reps,
    {
      indices_star <- sample(nrow(x), nrow(x), replace = TRUE)
      args_star <- bp_aggr_prepare(x[indices_star, ], values, kind)
      est_star <- do.call(bp_aggr_est_matrix, args_star)
      sd_star <- sqrt(do.call(bp_aggr_var_matrix, args_star))
      (transformer$est(est_star) - est_t) /
        transformer$sd(est_star, sd_star)
    },
    future.seed = TRUE
  )
}
