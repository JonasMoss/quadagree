x <- rbind(
  c(0, 0, 0, 0, 6),
  c(0, 0, 0, 0, 6),
  c(0, 0, 0, 0, 6)
)


#x <- fleissci::dat.fleiss1971
values = seq(ncol(x))
r <- sum(x[1, ])

xtx <- apply(x, 1, \(row) sum(values ^ 2 * row))
xt1 <- apply(x, 1, \(row) sum(values * row))
xt12 <- xt1^2

sxtx <- mean(xtx)
sxt1 <- mean(xt1)

dis <- 2/(r - 1) * (sxtx - sxt1^2 / r - var(xt1))

1 - 2*dis/r^2
