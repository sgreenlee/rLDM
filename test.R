y <- rep(c(-1,1), each=50)
y <- as.factor(y)
X <- matrix(rnorm(200), ncol=2)
X[51:75,] <- X[51:75,] + matrix(rep(c(2.5,-.5), each=25), ncol=2)
X[76:100,] <- X[76:100,] + matrix(rep(c(-1.5, 1.5), each=25), ncol=2)

