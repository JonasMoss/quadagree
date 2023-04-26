fleiss_prepare <- \(x, type) {
  calc <- list(xx = x, type = type)
}

fleiss_fun <- \(calc) {
  list(est = fleiss(calc$xx), var = avar(calc$xx, calc$type, TRUE))
}

conger_prepare <- \(x, type) {
  calc <- list(xx = x, type = type)
}

conger_fun <- \(calc) {
  list(est = conger(calc$xx), var = avar(calc$xx, calc$type, FALSE))
}

bp_prepare <- \(x, values, kind, type) {
  x <- as.matrix(x)
  if (is.null(values)) values <- unique(c(x))
  c1 <- bp_get_c1(values, kind)
  calc <- list(xx = x, c1 = c1, type = type)
}

pb_fun <- \(calc) {
  list(est = bp(calc$xx, calc$c1), var = avar(calc$xx, calc$type, calc$c1))
}
