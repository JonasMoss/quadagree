test_that("population and x-based measures coincide", {

  x = dat.zapf2016
  n = nrow(x)
  sigma = stats::cov(x) * (n - 1) / n
  mu = colMeans(x)

  expect_equal(
    fleiss(x),
    fleiss_pop(mu, sigma)
  )
  expect_equal(
    conger(x),
    conger_pop(mu, sigma)
  )
})
