#' @rdname fleiss
#' @export
fleiss_aggr <- \(x, values = seq_len(ncol(x))) {
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

#' @export
#' @rdname quadagree
bp_aggr <- function(x, values = seq_len(ncol(x)), kind = 1) {
  stopifnot(kind == 1 | kind == 2)
  stopifnot(ncol(x) == length(values))
  calc <- bp_aggr_prepare(x, values, kind)
  bp_aggr_est(calc)
}
