tune_ldm <- function(x, ...) UseMethod("tune_ldm")

tune_ldm.default <- function(X, y, nfolds=5) {
  
  ### Only works with radial kernels A.T.M. ####
  
  gamma<- c(2^-6, 2^-4, 2^-2, 2^-1, 1)
  lambda_1 <- seq(0,1, by=.1)
  lambda_2 <- seq(0,1, by=.1)
  cost <- 1:10

  folds <- sample(1:nfolds, nrow(X), replace=TRUE)
  
  params <- expand.grid(gamma, lambda_1, lambda_2, cost)
  n <- nrow(param_matrix)
  
  cv.errs <- vector(mode="numeric", length=n)
  min.err <- 1e7
  min.err.index <- 0
  for (i in 1:n) {
    cv.err <- 0
    for (k in 1:nfolds) {
      model <- ldm(X[folds!=k, ], y[folds!=k], 
                   gamma=params[i, 1], lambda_1=params[i, 2], lambda_2=params[i, 3], cost=params[i,4])
      pred <- predict(model, X[folds==k])
      cv.err <- cv.err + mean(pred==y[folds==k])
    }
    cv.err <- cv.err / nfolds
    cv.errs[i] <- cv.err
    if (cv.err < min.err) {
      min.err <- cv.err
      min.err.index <- i
    }
  }
  return(params[min.err.index, ])
}

tune_ldm.formula <- function(formula, data, ...) {
  
}