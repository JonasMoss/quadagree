test_that("population and x-based measures coincide", {
  x <- dat.zapf2016
  n <- nrow(x)
  values <- sort(unique(c(as.matrix(x))))
  sigma <- stats::cov(x) * (n - 1) / n
  mu <- colMeans(x)

  expect_equal(
    fleiss(x),
    fleiss_pop(mu, sigma)
  )
  expect_equal(
    conger(x),
    conger_pop(mu, sigma)
  )
  expect_equal(
    cohen(x),
    cohen_pop(mu, sigma)
  )
  expect_equal(
    bp(x),
    bp_pop(mu, sigma, values = values)
  )
})

test_that("agreement with irrCAC where relevant", {
  x <- dat.zapf2016
  est_irr <- (irrCAC::bp.coeff.raw(x, weights = "quadratic"))$est[4]
  est_quad <- bp(x, type = 1)
  expect_equal(as.numeric(est_irr), est_quad)
})
