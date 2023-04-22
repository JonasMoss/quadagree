#' Confidence intervals for the quadratically weighted Fleiss' kappa and
#'    Conger's (Cohen's) kappa
#'
#' Calculate confidence intervals quadratically weighted multi-rater kappas
#'    with potentially missing data. `congerci` and `cohenci` are aliases.
#'    See [fleiss()] and [cohen()] for estimation without inference.
#'
#' The confidence intervals are based on the formulas of Moss and van Oest
#'    (work in progress) along with standard asymptotic theory
#'    (Magnus, Neudecker, 2019) and the missing data theory of
#'    van Praag et al. (1985).
#'
#' The methods handle missing data using pairwise available information, i.e.,
#'    the option `use = "pairwise.complete.obs"` in [stats::cov()] along with
#'    the asymptotic theory of van Praag et al. (1985). The bootstrap option
#'    uses the studentized bootstrap (Efron, B. 1987), which is second order
#'    correct. Both functions makes use of [`future.apply`] when bootstrapping.
#'
#' The `type` variables defaults to `adf`, asymptotically distribution-free,
#'    which is consistent when the fourth moment is finite. The `normal` option
#'     assumes normality, and is not consistent for models with excess
#'    kurtosis unequal to `0`. The `elliptical` option assumes an
#'    elliptical or pseudo-elliptical distribution of the data. The resulting
#'    confidence intervals are corrected variants of the normal theory
#'    intervals with a kurtosis correction (Yuan & Bentler 2002). The
#'    common kurtosis parameter is calculated using the unbiased sample
#'    kurtosis (Joanes, 1998).
#'
#' Conger's (1980) kappa is a multi-rater generalization of Cohen's kappa.
#'    All functions in this package work for multiple raters, so functions
#'    starting with `cohen` or `conger` are aliases. The quadratically
#'    weighted Cohen's kappa is also known as Lin's concordance coefficient.
#'
#' The only difference between Cohen's kappa and Fleiss' kappa lies on how they
#'    measure disagreement due to chance. Here Fleiss' marginalizes the rating
#'    distribution across raters, essentially assuming there is no difference in
#'    the rating distribution across raters, while Cohen's kappa does not.
#'    There is a large literature comparing Fleiss' kappa to Cohen's kappa, and
#'    there is no consensus on which to prefer.
#'
#' @references
#'
#' Efron, B. (1987). Better Bootstrap Confidence Intervals. Journal of the
#' American Statistical Association, 82(397), 171-185.
#' https://doi.org/10.2307/2289144
#'
#' Van Praag, B. M. S., Dijkstra, T. K., & Van Velzen, J. (1985).
#' Least-squares theory based on general distributional assumptions with
#' an application to the incomplete observations problem.
#' Psychometrika, 50(1), 25–36. https://doi.org/10.1007/BF02294145
#'
#' Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample skewness
#' and kurtosis. Journal of the Royal Statistical Society: Series D
#' (The Statistician), 47(1), 183-189. https://doi.org/10.1111/1467-9884.00122
#'
#' Cohen, J. (1968). Weighted kappa: Nominal scale agreement with provision for
#' scaled disagreement or partial credit. Psychological Bulletin, 70(4),
#' 213–220. https://doi.org/10.1037/h0026256
#'
#' Fleiss, J. L. (1975). Measuring agreement between two judges on the presence
#' or absence of a trait. Biometrics, 31(3), 651–659.
#' https://www.ncbi.nlm.nih.gov/pubmed/1174623
#'
#' Conger, A. J. (1980). Integration and generalization of kappas for multiple
#' raters. Psychological Bulletin, 88(2), 322–328.
#' https://doi.org/10.1037/0033-2909.88.2.322
#'
#' Lin, L. I. (1989). A concordance correlation coefficient to evaluate
#' reproducibility. Biometrics, 45(1), 255–268.
#' https://www.ncbi.nlm.nih.gov/pubmed/2720055
#'
#' Moss, van Oest (work in progress). Inference for quadratically weighted
#' multi-rater kappas with missing raters.
#'
#' Moss (work in progress). On the Brennan–Prediger coefficients.
#'
#' Magnus, J. R., & Neudecker, H. (2019). Matrix Differential Calculus with
#' Applications in Statistics and Econometrics. John Wiley & Sons.
#' https://doi.org/10.1002/9781119541219
#'
#' @export
#' @param x Input data data can be converted to a matrix using `as.matrix`.
#' @param values to attach to each column on the Fleiss form data.
#'    Defaults to `1:C`, where `C` is the number of categories. Only used
#'    in `fleiss_aggr` and `bp_aggr`.
#' @param type Type of confidence interval. Either `adf`, `elliptical`, or
#'   `normal`. Ignored in `fleiss_aggrci`.
#' @param kind The kind of Brennan-Prediger coefficient used, `1` for the
#'   classical kind and `2` for the kind introduced in Moss (2023). Only
#'   relevant for `bp_aggr` and `bp`.
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
fleiss_aggr_ci <- function(x,
                           values = seq(ncol(x)),
                           transform = "none",
                           conf_level = 0.95,
                           alternative = c("two.sided", "greater", "less"),
                           bootstrap = FALSE,
                           n_reps = 1000) {
  call <- match.call()
  r <- sum(x[1, ])
  stopifnot(ncol(x) == length(values))
  alternative <- match.arg(alternative)
  transformer <- get_transformer(transform)

  quants <- limits(alternative, conf_level)
  x <- as.matrix(x)
  est <- fleiss_aggr(x, values)
  sd <- sqrt(fleiss_aggr_var(x, values))

  ci <- if (!bootstrap) {
    ci_asymptotic(est, sd, nrow(x), transformer, quants)
  } else {
    ci_boot_aggr(
      x,
      values,
      est,
      sd,
      transformer,
      quants,
      n_reps
    )
  }

  names(ci) <- quants
  attr(ci, "conf_level") <- conf_level
  attr(ci, "alternative") <- alternative
  attr(ci, "type") <- "Aggregated Fleiss kappa"
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

#' @export
#' @rdname fleissci
bp_aggr <- \(x, values = seq(ncol(x)), kind = 1) {
  stopifnot(kind == 1 | kind == 2)
  args <- bp_aggr_prepare(x, values, kind)
  c(
    est = do.call(bp_aggr_est_matrix, args),
    var = do.call(bp_aggr_var_matrix, args)
  )
}

#' @export
#' @rdname fleissci
cohenci <- congerci

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
    print(
      c(
        kappa = at("estimate"),
        sd = at("sd")
      ),
      digits = digits
    )
  }
  invisible(x)
}
