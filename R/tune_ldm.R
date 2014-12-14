tune_ldm <- function(x, ...) UseMethod("tune_ldm")

tune_ldm.default <- function(X, y, kernel='radial', scale=T, ranges=NULL, nfolds=5) {
  
  ### Only works with radial kernels A.T.M. ####
  if (length(scale)==1) scale <- rep(scale, ncol(X))
  else if (length(scale) != ncol(X)) {
    warning("Data not scaled: length of 'scale' argument must be either 1 or the number of columns of X")
    scale <- FALSE
  }
  
  x.scale <- NULL
  if (any(scale)) {
    tmp <- scale(X[,scale])
    X[,scale] <- tmp
    x.scale <- attributes(tmp)[c("scaled:center", "scaled:scale")]
  }
  
  if(isNaN(ranges)){
  ranges <- list(
    gamma = 2 ^ seq(-2:2),
    lambda_1 = 2 ^ seq(-8, -2, by=1),
    lambda_2 = 2 ^ seq(-8, -2, by=1),
    cost = c(10,50,100))
  }
  
  folds <- sample(1:nfolds, nrow(X), replace=TRUE)
  
  pmat <- expand.grid(ranges)
  n <- nrow(pmat)
  
  cv.errs <- vector(mode="numeric", length=n)
  min.err <- 1
  j <- 0
  for (i in 1:n) {
    cv.err <- 0
    for (k in 1:nfolds) {
      model <- do.call(ldm, c(list(X=X[folds!=k, ], y=y[folds!=k], kernel=kernel, scale=F), as.list(pmat[i,])))
      pred <- predict(model, X[folds==k,])
      cv.err <- cv.err +  (1 - mean(pred==factor(y[folds==k])))
    }
    cv.err <- cv.err / nfolds
    cv.errs[i] <- cv.err
    if (cv.err < min.err) {
      min.err <- cv.err
      j <- i
    }
  }
  best.mod <- do.call(ldm, c(list(X=X, y=y, scale=scale, kernel=kernel), as.list(pmat[j, ])))
  return(c(list(best.model=best.mod, cv.err=cv.errs[j]), as.list(pmat[j, ]))) 
}

tune_ldm.formula <- function(formula, data, ...) {
    mf <- model.frame(formula, data)
    Terms <- attr(mf, "terms")
    attr(Terms, "intercept") <- 0
    
    X <- model.matrix(Terms, mf)
    y <- model.extract(mf, "response")
    
    ret <- tune_ldm.default(X, y, ...)
    ret$call <- call
    ret$formula <- formula
    class(ret$best.model) <- c("ldm.formula", "ldm")
    return(ret)
}