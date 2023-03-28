### ===========================================================================
###   Testing weighted least squares.
### ===========================================================================

data <- fleissci::dat.zapf2016
sampler <- \(n) dplyr::slice_sample(data, n = n)

wls <- \(z) {
  yy <- construct_y(z)
  lmmod <- lm(yy$y1 ~ yy$y2)
  wlmmod <- lm(yy$y1 ~ yy$y2, w = yy$w)
  c(coef(lmmod)[2],
    coef(wlmmod)[2],
    summary(lmmod)$coefficients[2, 2],
    summary(wlmmod)$coefficients[2, 2])
}

results <- t(replicate(5000, {
  z <- sampler(20)
  wls(z)
}))

mean(results[, 1])
mean(results[, 2])
fleissci::fleiss(data)


sd(results[, 1])
mean(results[, 3])

sd(results[, 2])
mean(results[, 4])

hist(results[, 2], breaks = 100)
hist(results[, 1], breaks = 100)
