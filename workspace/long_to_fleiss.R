#' Convert data on raw form to Fleiss form.
#'
#' @param x Data on raw form.
#' @param n_cat Number of categories.
#' @return Data on Fleiss form.
raw_to_fleiss <- \(x, n_cat) {
  n <- nrow(x)
  y <- matrix(0, nrow = n, ncol = n_cat)
  for (i in seq(n)) {
    for (j in x[i, ]) {
      y[i, j] = 1 + y[i, j]
    }
  }
  y
}

#' Transform data Fleiss form to "fake data" on wide form.
#'
#' Data on Fleiss form doesn't keep which judge made which rating. Instead, it
#'    aggregates them into counts. This function disaggregates the data into
#'    one of the possible sources. There are many, but that doesn't
#'    matter for the computation of Fleiss' kappa.
#'
#' @param x Data on the Fleiss form.
#' @keywords internal
#' @return Fake wide form data.
fleiss_to_wide <- function(x) {
  x <- as.matrix(x)
  y <- matrix(nrow = nrow(x), ncol = sum(x[1, ]))
  for (row in seq(nrow(y))) {
    current <- 1
    for (col in seq(ncol(y))) {
      while (TRUE) {
        if (x[row, current] != 0) {
          y[row, col] <- current
          x[row, current] <- x[row, current] - 1
          break
        } else {
          current <- current + 1
        }
      }
    }
  }
  y
}
