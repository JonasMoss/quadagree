A = matrix(c(2,1,1,2), 2, 2)
# positive definite
eig = eigen(A)
d = eig$values
q = eig$vectors
r <- seq(0, 2*pi, 0.01)
zz = rbind(1 / sqrt(d[1])*sin(r), 1/sqrt(d[2])*cos(r))
plot(t(zz))
xx <- q %*% zz
plot(t(xx), type = "l", lwd = 2)
sapply(seq_along(r), \(i) t(t(xx)[i, ]) %*% A %*% t(xx)[i, ])

A = matrix(c(2,2,2,2), 2, 2)
# positive semi-definite
eig = eigen(A)
d = eig$values
q = eig$vectors
r <- seq(-10, 10, 0.01)
zz = rbind(1/sqrt(d[1]), r)
plot(t(zz))
xx <- q %*% zz
plot(t(xx), type = "l", lwd = 2)
sapply(seq_along(r), \(i) t(t(xx)[i, ]) %*% A %*% t(xx)[i, ])





A = matrix(c(1,2,2,2), 2, 2)
# not pos/neg definite
eig = eigen(A)
d = eig$values
q = eig$vectors

r <- c(seq(-pi/2+0.1, pi/2-0.1, 0.01), seq(pi-pi/2+0.1, pi + pi/2-0.1, 0.01))
zz = rbind(1/sqrt(d[1]) / cos(r), 1 / sqrt(-d[2]) * tan(r))
plot(t(zz))

xx <- q %*% zz
plot(t(xx), type = "l", lwd = 2)
sapply(seq_along(r), \(i) t(t(xx)[i, ]) %*% A %*% t(xx)[i, ])



## Rotation

r = pi/3
A =matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2, 2)
A =matrix(c(1, 0, 1, 1), 2, 2)



B =matrix(c(0, 0, 1, 0), 2, 2)
A = matrix(c(2,1,1,2), 2, 2)
B =matrix(c(-1, 0, 1, 0), 2, 2)
B =matrix(c(0, -1, 1, 0), 2, 2)
eigen(B)

B = -matrix(c(1, 0, 0, -1), 2, 2)
diag(sqrt(c(1 + 0i,-1 + 0i)))


N = matrix(c(0, 0, 1, 0), 2, 2)
S = 2 * diag(2)

A = N + S
Asym = 0.5 * ((N + S) + t(N+S))
Askew = 0.5 * ((N + S) - t(N+S))

Asym %*% Askew - Askew %*% Asym
1/4 * (Asym %*% Askew + Askew %*% Asym)


B = matrix(c(1, 2, 3, 2, 3, 1, 1, 1, 2),3,3)
B = matrix(runif(9),3,3)

B_sym = t(lower.tri(B, diag = TRUE) * B) + lower.tri(B, diag = FALSE) * B
N = B - B_sym

A = matrix(c(1, 0, 1, 1, 1, 0, 0, 1, 1), 3, 3)
A_sym = t(upper.tri(A, diag = TRUE) * A) + upper.tri(A, diag = FALSE) * A
A - A_sym
