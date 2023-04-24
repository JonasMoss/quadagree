#' @export
#' @rdname quadagree
fleissci_aggr <- function(x,
                          values = seq_len(ncol(x)),
                          transform = "none",
                          conf_level = 0.95,
                          alternative = c("two.sided", "greater", "less"),
                          bootstrap = FALSE,
                          n_reps = 1000) {
  call <- match.call()
  stopifnot(ncol(x) == length(values))
  calc <- fleiss_aggr_prepare(x, values)
  est_fun <- fleiss_aggr_est
  var_fun <- fleiss_aggr_var
  quadagree_aggr_(
    calc = calc,
    transform = transform,
    conf_level = conf_level,
    alternative = alternative,
    bootstrap = bootstrap,
    n_reps = n_reps,
    est_fun = est_fun,
    var_fun = var_fun,
    call = quote(call)
  )
}

#' @export
#' @rdname quadagree
bpci_aggr <- function(x,
                      values = seq_len(ncol(x)),
                      kind = 1,
                      transform = "none",
                      conf_level = 0.95,
                      alternative = c("two.sided", "greater", "less"),
                      bootstrap = FALSE,
                      n_reps = 1000) {
  stopifnot(kind == 1 | kind == 2)
  stopifnot(ncol(x) == length(values))
  call <- match.call()
  quadagree_aggr_(
    calc = bp_aggr_prepare(x, values, kind),
    transform = transform,
    conf_level = conf_level,
    alternative = alternative,
    bootstrap = bootstrap,
    n_reps = n_reps,
    est_fun = bp_aggr_est,
    var_fun = bp_aggr_var,
    call = quote(call)
  )
}

#' @keywords internal
quadagree_aggr_ <- function(calc,
                            transform = "none",
                            conf_level = 0.95,
                            alternative = c("two.sided", "greater", "less"),
                            bootstrap = FALSE,
                            n_reps = 1000,
                            est_fun,
                            var_fun,
                            call) {
  alternative <- match.arg(alternative)
  transformer <- get_transformer(transform)

  quants <- limits(alternative, conf_level)
  est <- est_fun(calc)
  sd <- sqrt(var_fun(calc))

  ci <- if (!bootstrap) {
    ci_asymptotic(est, sd, calc$n, transformer, quants)
  } else {
    ci_boot_aggr(
      est,
      sd,
      calc,
      transformer,
      quants,
      n_reps,
      est_fun = est_fun,
      var_fun = var_fun
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
