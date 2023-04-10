x <- fleissci::dat.fleiss1971
fleiss_aggr_ci(x, bootstrap = TRUE)
fleiss_aggr_ci(x)

fleissci(fleiss_to_wide(x))


x <- agreeable::simulate_jsm(100, c(0.6, 0.7, 0.8), model = "fleiss",
                        true_dist = c(0.1,0.1,0.5,0.3))
skill <- attr(x, "skill")

results <- replicate(10000, {
  x <- agreeable::simulate_jsm(100, c(0.6, 0.7, 0.8), model = "fleiss",
                               true_dist = c(0.1,0.1,0.5,0.3))
  x <- raw_to_fleiss(x, n_cat = 4)
  ci <- fleiss_aggr_ci(x, bootstrap = TRUE)
  ci[1] < skill & ci[2] > skill
})

mean(results)
