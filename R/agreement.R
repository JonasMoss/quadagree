#' Calculate the quadratic Fleiss' kappa and Conger's (Cohen's) kappa
#'
#' Estimate a quadratically weighted kappa with potentially missing data using
#'    least squares, or calculate its population value. `conger` and `cohen` are
#'    aliases. See [quadagree()] and [cohenci()] for confidence intervals.
#'
#' Conger's kappa is a multi-rater generalization of Cohen's kappa.
#'    All functions in this package work for multiple raters, so functions
#'    starting with `cohen` or `conger` are aliases. The quadratically weighted
#'    Cohen's kappa is also known as Lin's (1989) concordance coefficient.
#'
#' The only difference between Cohen's kappa and Fleiss' kappa lies on how they
#'    measure disagreement due to chance. Here Fleiss' marginalizes the rating
#'    distribution across raters, essentially assuming there is no difference in
#'    the rating distribution across raters, while Cohen's kappa does not.
#'    There is a large literature comparing Fleiss' kappa to Cohen's kappa,
#'    and there is no consensus on which to prefer.
#'
#' The functions `cohen_pop`, `conger_pop`, and `fleiss_pop` calculate the
#'    population values of the kappas using the mean and covariance of the
#'    rating distributions. The derivation of the covariance and mean
#'    formulation can be found in (Moss and van Oest, work in progress.)
#'
#' The functions `cohen`, `conger`, and `fleiss` estimate the kappas using least
#'    squares. This method uses the biased sample covariance estimator, which
#'    agrees with the literature on agreement coefficients in the case of two
#'    raters. In the case of missing data, the pairwise available data is used,
#'    employing the option `use = "pairwise.complete.obs"` in [stats::cov()].
#'
#'
#' @param x Input data on wide form. Columns are raters, items are ratings.
#'    Missing values are allowed, and should be encoded as `NA`s.
#' @param mu Population vector of means.
#' @param sigma Population covariance matrix.
#' @param values to attach to each column on the Fleiss form data.
#'    Defaults to `1:C`, where `C` is the number of categories. Only used
#'    in `fleiss_aggr`.
#' @return Sample / population quadratically weighted Fleiss' or Conger's kappa.
#' @examples
#' x <- irrCAC::cac.raw4raters[2:9, ]
#' fleiss(x)
#' # [1] 0.6666667
#' irrCAC::fleiss.kappa.raw(x, weights = "quadratic")$est$coeff.val
#' # [1] 0.66667 (irrCAC prematurely rounds the coefficient)
#'
#' x <- irrCAC::cac.raw4raters[2:9, ]
#' conger(x)
#' # [1] 0.6719243
#' irrCAC::conger.kappa.raw(x, weights = "quadratic")$est$coeff.val
#' # [1] 0.67192 (irrCAC prematurely rounds the coefficient)
#' @references
#' Cohen, J. (1968). Weighted kappa: Nominal scale agreement with provision for
#'  scaled disagreement or partial credit. Psychological Bulletin, 70(4),
#'  213–220. https://doi.org/10.1037/h0026256
#'
#' Fleiss, J. L. (1975). Measuring agreement between two judges on the presence
#'  or absence of a trait. Biometrics, 31(3), 651–659.
#' https://www.ncbi.nlm.nih.gov/pubmed/1174623
#'
#' Conger, A. J. (1980). Integration and generalization of kappas for
#' multiple raters. Psychological Bulletin, 88(2), 322–328.
#' https://doi.org/10.1037/0033-2909.88.2.322
#'
#' Lin, L. I. (1989). A concordance correlation coefficient to evaluate
#' reproducibility. Biometrics, 45(1), 255–268.
#' https://www.ncbi.nlm.nih.gov/pubmed/2720055
#'
#' Moss, van Oest (work in progress). Inference for quadratically weighted
#' multi-rater kappas with missing raters.
#'
#' @export

#' @rdname fleiss
#' @export
bp <- \(x, values = NULL, type = 1) {
  y <- as.matrix(x)

  if (is.null(values)) {
    values <- sort(unique(c(y)))
  }

  n <- nrow(y)

  sigma <- stats::cov(y, use = "pairwise.complete.obs") * (n - 1) / n

  if (any(is.na(sigma))) {
    stop("The data does not contain sufficient non-NAs.")
  }

  mu <- colMeans(y, na.rm = TRUE)
  bp_pop(mu, sigma, values, type)
}

#' @name fleiss
fleiss <- function(x) {
  n <- nrow(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) {
    stop("The data does not contain sufficient non-NAs.")
  }
  mu <- colMeans(x, na.rm = TRUE)
  fleiss_pop(mu, sigma)
}

#' @rdname fleiss
#' @export
conger <- function(x) {
  n <- nrow(x)
  sigma <- stats::cov(x, use = "pairwise.complete.obs") * (n - 1) / n
  if (any(is.na(sigma))) {
    stop("The data does not contain sufficient non-NAs.")
  }
  mu <- colMeans(x, na.rm = TRUE)
  conger_pop(mu, sigma)
}

#' @rdname fleiss
#' @export
cohen <- function(x) conger(x)

#' @rdname fleiss
#' @export
bp_pop <- function(mu, sigma, values, type = 1) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  mean_diff <- (mean(mu^2) - mean(mu)^2) * r
  mean_sum <- sum(sigma) / r
  c1 <- bp_aggr_get_c1(values, type)
  d <- 2 / (r - 1) * (trace + mean_diff - mean_sum)
  1 - d / c1
}

#' @rdname fleiss
#' @export
conger_pop <- function(mu, sigma) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  top <- sum(sigma) - trace
  bottom <- (r - 1) * trace + r^2 * (mean(mu^2) - mean(mu)^2)
  top / bottom
}

#' @rdname fleiss
#' @export
cohen_pop <- function(mu, sigma) conger_pop(mu, sigma)

#' @rdname fleiss
#' @export
fleiss_pop <- function(mu, sigma) {
  r <- ncol(sigma)
  trace <- sum(diag(sigma))
  top <- sum(sigma) - trace - r * (mean(mu^2) - mean(mu)^2)
  bottom <- (r - 1) * trace + (r - 1) * r * (mean(mu^2) - mean(mu)^2)
  top / bottom
}

#' @rdname fleiss
#' @export
fleiss_aggr <- \(x, values = seq(ncol(x))) {
  r <- sum(x[1, ])
  stopifnot(ncol(x) == length(values))

  y <- as.matrix(x)
  xtx <- tcrossprod(values^2, y)
  xt1 <- tcrossprod(values, y)

  extx <- mean(xtx)
  ext1 <- mean(xt1)
  ext2 <- mean(xt1^2)

  1 / (r - 1) * ((ext2 - ext1^2) / (extx - ext1^2 / r) - 1)
}
