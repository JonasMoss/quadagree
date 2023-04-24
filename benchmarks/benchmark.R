x = dat.fleiss1971
x = rbind(x)
microbenchmark::microbenchmark(
  irrCAC::bp.coeff.dist(x, weights = "quadratic"),
  bpci_aggr(x),
  irrCAC::fleiss.kappa.dist(x, weights = "quadratic"),
  fleissci_aggr(x)
)

irrCAC::fleiss.kappa.dist(x, weights = "quadratic")
attr(fleissci_aggr(x), "est")
attr(fleissci_aggr(x), "sd")/sqrt(nrow(x) - 1)

irrCAC::bp.coeff.dist(x, weights = "quadratic")
attr(bpci_aggr(x), "est")
attr(bpci_aggr(x), "sd")/sqrt(nrow(x) - 1)

x = dat.zapf2016
x = rbind(x, x, x, x, x, x, x, x, x)
microbenchmark::microbenchmark(
  irrCAC::conger.kappa.raw(x, weights = "quadratic"),
  congerci(x),
  irrCAC::fleiss.kappa.raw(x, weights = "quadratic"),
  fleissci(x)
)

irrCAC::fleiss.kappa.dist(x, weights = "quadratic")
attr(fleissci_aggr(x), "est")
attr(fleissci_aggr(x), "sd")/sqrt(nrow(x) - 1)

irrCAC::bp.coeff.dist(x, weights = "quadratic")
attr(bpci_aggr(x), "est")
attr(bpci_aggr(x), "sd")/sqrt(nrow(x) - 1)

