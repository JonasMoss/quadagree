data <- fleissci::dat.zapf2016

construct_y <- \(z, vars = NULL) {
  r <- ncol(z)
  n <- nrow(z)
  indices <- arrangements::combinations(r, 2)
  indices <- rbind(indices, cbind(indices[, 2], indices[, 1]))

  y1 <- y2 <- c()
  index1 <- index2 <- c()
  for (row in seq(nrow(indices))) {
    y1 <- c(y1, z[, indices[row, 1]])
    y2 <- c(y2, z[, indices[row, 2]])
    index1 <- c(index1, rep(indices[row, 1], n))
    index2 <- c(index2, rep(indices[row, 2], n))
  }

  uniques <- sort(unique(c(as.matrix(z))))
  if(is.null(vars)) {
    vars <- sapply(uniques, \(x) var(y1[y2 == x]))
    vars[vars == 0] = mean(vars[vars != 0])
    names(vars) <- uniques
  }

  covs <- outer(uniques, uniques, Vectorize(\(x1, x2) covar(x1, x2, z)))
  colnames(covs) <- rownames(covs) <- uniques

  ii <- rep(seq(nrow(z)), r * (r-1))

  w <- matrix(0, length(y1), length(y1))
  trans <- sapply(y2, \(x) which(uniques == x))
  for (i in seq_along(y1)) {
    for (j in seq_along(y1)) {
      if(i == j) {
        0
      } else {
        # if (ii[i] == ii[j]) {
        #   w[i,j] = covs[trans[i], trans[j]]
        # }
        w[i,j] = covs[trans[i], trans[j]]/n
      }
    }
  }

  diag(w) = vars[sapply(y1, \(x) which(names(vars) == x))]


  list(y1 = y1, y2 = y2, w = w)
}

z <- data
yy <- construct_y(z)
ii <- rep(seq(nrow(z)), r * (r-1))

summary(lm(yy$y1 ~ yy$y2 + ii))
summary(lm(yy$y1 ~ yy$y2))
#summary(lm(yy$y1 ~ yy$y2, w = yy$w))

MASS::lm.gls(yy$y1 ~ yy$y2, W = solve(yy$w))$coefficients


xx = cbind(1, yy$y2)
solve(t(xx) %*% solve(yy$w) %*% xx)[2, 2]


cor(yy$y1, yy$y2, method = "kendall")

table(cbind(yy$y1, yy$y2))


plot(table(yy$y1, yy$y2))


rr = arrangements::combinations(4, 2)
mean(apply(rr, 1, \(r) cor(c(data[, r[1]], data[, r[2]]), c(data[, r[2]], data[, r[1]]), method = "kendall")))
mean(apply(rr, 1, \(r) cor(c(data[, r[1]], data[, r[2]]), c(data[, r[2]], data[, r[1]]))))
cor(yy$y1, yy$y2, method = "kendall")
cor(yy$y1, yy$y2)


tab = table(yy$y1, yy$y2)
probs = tab/rbind(colSums(tab), colSums(tab),colSums(tab),colSums(tab),colSums(tab))
means = sapply(1:5, \(x) mean(as.matrix(data) == x))

sum(diag(tab) / rowSums(tab) * means)


sqrt(summary(lm(yy$y1 ~ as.factor(yy$y2)))$r.sq)
irrCAC::fleiss.kappa.raw(data)

L1pack::lad(yy$y1 ~ yy$y2)
L1pack::lad(yy$y1 ~ as.factor(yy$y2))

#L1pack::lad(c(data[, 1], data[, 2]) ~ c(data[, 2], data[, 1]))
#fitted <- mod$fitted.values
x <- c(data[, 2], data[, 1])
y <- c(data[, 1], data[, 2])
1 - mean(abs(x - y)) / mean(abs(x - median(x)))
irrCAC::fleiss.kappa.raw(data[, c(1, 2)], weights = "linear")



x_ <- data[, 1] + rnorm(50, 0, 0.001)
y_ <- data[, 2] + rnorm(50, 0, 0.001)
x <- c(x_, y_)
y <- c(y_, x_)
cor(x_, y_, method = "kendall")
cor(x, y, method = "kendall")


table(as.matrix(data))

mean(as.matrix(data) == 5)^2 +
  mean(as.matrix(data) == 4)^2 +
  mean(as.matrix(data) == 3)^2 +
  mean(as.matrix(data) == 2)^2 +
  mean(as.matrix(data) == 1)^2

irrCAC::fleiss.kappa.raw(data)

data1 = data == 5
irrCAC::fleiss.kappa.raw(data1)

yy_ <- construct_y(data1)
lm(yy_$y1 ~ yy_$y2)

### Test

pis = sapply(1:5, \(i) mean(as.matrix(data) == i))


data = data[, 1:2]
cors = sapply(1:5, \(i) {
  datai = data == i
  yy_ = construct_y(datai)
  cor(yy_$y1, yy_$y2)
})

vars = sapply(1:5, \(i) {
  datai = data == i
  yy_ = construct_y(datai)
  var(yy_$y1,)
})


sum(cors * vars)/sum(pis^2)
irrCAC::fleiss.kappa.raw(data)


fleissci::fleiss(datai[, c(1,2)])
irrCAC::fleiss.kappa.raw(datai[, c(1,2)])


y1 <- yy$y1
y2 <- yy$y2

eig <- eigen(w)
positives <- (eig$values > 0)
t1 <- c(eig$vectors %*% y1)[positives]
t2 <- c(eig$vectors %*% y2)[positives]
d_inv <- diag(1/eig$values[positives])

t(t2) %*% d_inv %*% t2

summary(lm(y1 ~ y2))

attr(fleissci::fleissci(data), "sd")/sqrt(50)

fleissci::fleissci(data[1:10, c(1,2,3)])
attr(fleissci::fleissci(data[1:10, c(1,2,3)]), "sd")/sqrt(10)




var_ = \(x) {
  result = var(x) * (length(x) - 1) / length(x)
  if(is.na(result)) 0 else result
}

p1 = colSums(data == 1)/sum(data == 1)
means = sapply(1:4, \(i) colMeans(data[,-i]))
vars = apply(means, 2, var_)

mean(sapply(1:4, \(i) var_(data[, i]) * p1[i]))

indices <- arrangements::combinations(r, 2)
x1 = data[, 1]
x2 = data[, 2]
x3 = data[, 3]
x4 = data[, 4]

a = 4
b = 5
cov(x1[x3 == a & x4 == b], x2[x3 == a & x4 == b])


## This is wrong - not needed.
a = 5
(var_(x1[x2 == a]) + var_(x1[x3 == a]) + var_(x1[x4 == a]) +
var_(x2[x1 == a]) + var_(x2[x3 == a]) + var_(x2[x4 == a]) +
var_(x3[x2 == a]) + var_(x3[x1 == a]) + var_(x3[x4 == a]) +
var_(x4[x2 == a]) + var_(x3[x3 == a]) + var_(x4[x1 == a])) / (4 * 3)

# for 2 and 1.

(data[, 1] == 1) * (data[, 1] == 1)
