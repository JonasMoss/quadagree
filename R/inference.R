#' Confidence intervals for quadratically weighted Fleiss' kappa and
#'    Conger's kappa.
#'
#' @export
#' @param x Input data data can be converted to a matrix using `as.matrix`.
#' @param type Type of confidence interval. Either `adf`, `elliptical`, or
#'   `normal`.
#' @param transform One of `"none"`, `"log"`, `"fisher"`, and `"arcsin`.
#'   Defaults to `"none"`.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of `"two.sided"` (default), `"greater"` or `"less"`.
#' @param conf_level Confidence level. Defaults to `0.95`.
#' @param bootstrap If `TRUE`, performs a studentized bootstrap with `n_reps`
#'   repetitions. Defaults to `FALSE`.
#' @param n_reps Number of bootstrap samples if `bootstrap = TRUE`. Ignored if
#'   `bootstrap = FALSE`. Defaults to `1000`.
#' @return A vector of class `fleissci` containing the confidence end points.
#'   The arguments of the function call are included as attributes.
#' @name fleissci

fleissci <- function(x,
                     type = c("adf", "elliptical", "normal"),
                     transform = "none",
                     conf_level = 0.95,
                     alternative = c("two.sided", "greater", "less"),
                     bootstrap = FALSE,
                     n_reps = 1000) {
  call <- match.call()
  args <- sapply(names(formals()), str2lang)
  do.call(what = fleissci_, c(args, call = quote(call), fleiss = TRUE))
}

#' @export
#' @rdname fleissci
congerci <- function(x,
                     type = c("adf", "elliptical", "normal"),
                     transform = "none",
                     conf_level = 0.95,
                     alternative = c("two.sided", "greater", "less"),
                     bootstrap = FALSE,
                     n_reps = 1000) {
  call <- match.call()
  args <- sapply(names(formals()), str2lang)
  do.call(what = fleissci_, c(args, call = quote(call), fleiss = FALSE))
}

#' @keywords internal
fleissci_ <- function(x,
                      type = c("adf", "elliptical", "normal"),
                      transform,
                      conf_level,
                      alternative = c("two.sided", "greater", "less"),
                      bootstrap,
                      n_reps,
                      call,
                      fleiss) {
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  transformer <- get_transformer(transform)

  quants <- limits(alternative, conf_level)
  x <- as.matrix(x)

  est <- if (!fleiss) conger(x) else fleiss(x)
  sd <- sqrt(avar(x, type, fleiss))

  ci <- if (!bootstrap) {
    ci_asymptotic(est, sd, nrow(x), transformer, quants)
  } else {
    ci_boot(
      x,
      est,
      sd,
      type,
      transformer,
      quants,
      n_reps,
      fleiss = fleiss
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
  class(ci) <- "fleissci"
  ci[2] <- min(ci[2], 1)
  ci
}

#' @export
print.fleissci <- function(x, digits = getOption("digits"), ...) {
  at <- \(y) attr(x, y)
  cat("Call: ", paste(deparse(at("call")),
    sep = "\n",
    collapse = "\n"
  ), "\n\n", sep = "")

  if (!is.null(x)) {
    cat(format(100 * at("conf_level")),
      "% confidence interval (n = ", at("n"), ").\n",
      sep = ""
    )
    print(x[1:2], digits = digits)
    cat("\n")
  }

  if (!is.null(at("estimate"))) {
    cat("Sample estimates.\n")
    print(c(
      kappa = at("estimate"),
      sd = at("sd")
    ),
    digits = digits
    )
  }
  invisible(x)
}
