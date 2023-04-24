true_dist <- c(0.8, 0.05, 0.05, 0.05, 0.05)
skills <- c(0.9, 0.9, 0.8, 0.9, 0.8, 0.7)

x <- agreeable::simulate_jsm(1000000, skills,
  model = "fleiss",
  true_dist = true_dist
)

avar(x, "adf", TRUE)


results <- replicate(200000, {
  x <- agreeable::simulate_jsm(
    2000,
    skills,
    model = "fleiss",
    true_dist = true_dist
  )
  fleiss(x)
})

var(results) * 2000
