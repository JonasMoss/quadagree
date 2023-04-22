true_value <- 0.6366667

xx <- agreeable::simulate_jsm(
  n = 100000,
  s = c(0.7, 0.5, 0.9),
  model = "fleiss",
  true_dist = c(0.2, 0.2, 0.3, 0.2, 0.1)
)
x <- raw_to_fleiss(xx, n_cat = 5)

fleiss_inference(x)


n_reps <- 10000
n <- 100
results <- replicate(n_reps, {
  xx <- agreeable::simulate_jsm(
    n = n,
    s = c(0.7, 0.5, 0.9),
    model = "fleiss",
    true_dist = c(0.2, 0.2, 0.3, 0.2, 0.1)
  )
  quadagree::fleiss(xx)
})

var(results) * n
