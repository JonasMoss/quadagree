#' @keywords internal
quadagree_ <- function(x,
                       type = c("adf", "elliptical", "normal"),
                       transform,
                       conf_level,
                       alternative = c("two.sided", "greater", "less"),
                       bootstrap,
                       n_reps,
                       call,
                       fun) {
  type <- match.arg(type)
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
    ci_boot(
      calc,
      fun,
      type,
      transformer,
      quants,
      n_reps
    )
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
quadagree_aggr_ <- function(calc,
                            transform = "none",
                            conf_level = 0.95,
                            alternative = c("two.sided", "greater", "less"),
                            bootstrap = FALSE,
                            n_reps = 1000,
                            fun,
                            call) {
  alternative <- match.arg(alternative)
  transformer <- get_transformer(transform)
  quants <- limits(alternative, conf_level)
  est_var <- fun(calc)
  est <- est_var$est
  sd <- sqrt(est_var$var)


  ci <- if (!bootstrap) {
    ci_asymptotic(est, sd, calc$n, transformer, quants)
  } else {
    ci_boot(
      calc,
      fun,
      type = NULL,
      transformer,
      quants,
      n_reps
    )
  }

  names(ci) <- quants
  attr(ci, "conf_level") <- conf_level
  attr(ci, "alternative") <- alternative
  attr(ci, "n") <- calc$n
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
