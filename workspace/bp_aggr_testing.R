### ============================================================================
###   Testing inference for aggregated Brennan-Prediger.
### ============================================================================

# Testing data
x <- fleissci::dat.fleiss1971

bp_aggr(x)
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
      y[i, j] <- 1 + y[i, j]
    }
  }
  y
}




microbenchmark::microbenchmark(
  bp(x, values = seq(ncol(x)), 1),
  irrCAC::bp.coeff.dist(x)
)



true_dist <- c(0.1, 0.1, 0.5, 0.3)
skills <- c(0.9, 0.9, 0.8)

x <- agreeable::simulate_jsm(100000, skills,
  model = "fleiss",
  true_dist = true_dist
)
x <- raw_to_fleiss(x, n_cat = length(true_dist))

results <- replicate(20000, {
  x <- agreeable::simulate_jsm(200, skills,
    model = "fleiss",
    true_dist = true_dist
  )
  x <- raw_to_fleiss(x, n_cat = length(true_dist))
  bp1(x, 1:4)
})

var(results) * 200
