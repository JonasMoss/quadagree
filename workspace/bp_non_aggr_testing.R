x <- dat.zapf2016[sample(50, 40), ]
n <- nrow(x)
r <- ncol(x)
c1 <- bp_aggr_get_c1(1:5, 1)

(irrCAC::bp.coeff.raw(x, weights = "quadratic")$est[5] * sqrt(n-1))^2
avar_bp(x, "adf", c1)


#(irrCAC::fleiss.kappa.raw(x, weights = "quadratic")$est[5] * sqrt(n-1))^2
#avar(x, "adf", TRUE)

s = c(0.9, 0.7, 0.6, 0.8, 0.1)
true_dist = c(0.3, 0.3, 0.3, 0.1)
c1 <- bp_aggr_get_c1(1:4, 1)
mod = agreeable::simulate_jsm(100000, s, model = "bp", true_dist = true_dist)
true = attr(mod, "skill")
c1 <- bp_aggr_get_c1(1:4, 1)
avar_bp(mod, "adf", c1)

res = replicate(10000, {
  x = agreeable::simulate_jsm(100, s, model = "bp", true_dist = true_dist)
  bp(x, values = 1:4, kind = 1)
})

var(res) * 100
