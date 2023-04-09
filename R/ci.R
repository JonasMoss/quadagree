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
                    fleiss = FALSE) {
  boots <- studentized_boots(n_reps, x, type, transformer)
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
                              fleiss = FALSE) {
  fun <- if (fleiss) fleiss else conger
  est <- fun(x)
  future.apply::future_replicate(n_reps,
    {
      indices_star <- sample(nrow(x), nrow(x), replace = TRUE)
      est_star <- fun(x)
      sd_star <- sqrt(avar(x[indices_star, ], type, fleiss))
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


#' Calculate limits of a confidence interval.
#'
#' @param alternative Alternative choosen.
#' @param conf_level Confidence level.
#' @keywords internal
limits <- function(alternative, conf_level) {
  half <- (1 - conf_level) / 2
  if (alternative == "two.sided") {
    return(c(half, 1 - half))
  }
  if (alternative == "greater") {
    return(c(2 * half, 1))
  }
  if (alternative == "less") {
    return(c(0, conf_level))
  }
}
