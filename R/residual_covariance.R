resid_fun <- \(z) {
  yy <- construct_y(z)
  wlmmod <- lm(yy$y1 ~ yy$y2, w = yy$w)
  resid(wlmmod)
}

z <- sampler(n)
yy <- construct_y(z)
indices <- rep(seq(n), r*(r-1))

resid_fun(z)


### Covariance

covar <- \(x1, x2, z) {
  r = ncol(z)

  covs = \(rr, x1, x2) {
    indices = (z[, rr[3]] == x1) * (z[, rr[4]] == x2) == 1
    cov(z[indices, rr[1]], z[indices, rr[2]])
  }

  expectations = \(rr, x) {
    indices = (z[, rr[2]] == x)
    if(sum(indices) == 0) NA else mean(z[indices, rr[1]])
  }


  rrs <- arrangements::permutations(r, 4)
  cov_mean <- mean(apply(rrs, 1, \(rr) covs(rr, x1, x2)))
  cov_mean <- if (is.na(cov_mean)) 0 else cov_mean

  rrs <- arrangements::permutations(r, 2)

  exp1 <- apply(rrs, 1, \(rr) {
    expectations(rr, x1)
  })

  exp2 <- apply(rrs, 1, \(rr) {
    expectations(rr, x2)
  })

  mean_cov <- cov(exp1, exp2, use = "pairwise.complete.obs")

  result <- cov_mean + mean_cov
  if(is.na(result)) 0 else result
}

outer(1:5, 1:5, Vectorize(\(x1, x2) covar(x1, x2, z)))
