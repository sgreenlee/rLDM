y <- rep(c(-1,1), each=50)
y <- as.factor(y)
X <- matrix(rnorm(200), ncol=2)
X[51:75,] <- X[51:75,] + matrix(rep(c(2.5,-.5), each=25), ncol=2)
X[76:100,] <- X[76:100,] + matrix(rep(c(-1.5, 1.5), each=25), ncol=2)
plot(X[51:100,], col=c('red'), pch=19)
points(X[1:50,], col='blue', pch=19)
setwd('~/projects/rldm')
debug(ldm)
fit.ldm <- ldm(y, X)
fit.svm <- svm(X, y, kernel='radial', cost=10, gamma=1)

table(y,fit.ldm$fitted) # training set performance for LDM
table(y,fit.svm$fitted) # training set performance for SVM

table(fit.ldm$fitted, fit.svm$fitted) # LDM fit vs SVM fit -- identical


