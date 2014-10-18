
ldm <- function(x, ...) UseMethod("ldm")

ldm.default <- function(y, X, kernel = 'radial', cost = 10, gamma=1, lambda_1=1/4, lambda_2=1/4)
{
  kernel_type <- pmatch(kernel, c('linear', 'polynomial', 'radial', 'sigmoid')) - 1
  
  if (!is.factor(y)) stop("y must be a factor")
  if (length(levels(y)) != 2) stop("y must be a two-level factor: multi-class not implemented yet")
  else y <- c(-1,1)[unclass(y)]
  
  X <- scale(X)
  
  ret <- ldm_train(y, X, gamma, lambda_1, lambda_2, cost)
  
  ret$cost <- cost
  ret$gamma <- gamma
  ret$lambda_1 <- lambda_1
  ret$lambda_2 <- lambda_2
  ret$call <- match.call()
  ret$fitted <- sign(ret$K %*% ret$alpha)
  class(ret) <- "ldm"
  
  ret
}
  
ldm.formula <- function(formula, data = NULL, ... )
{
  call <- match.call()
  if(is.matrix(eval.parent(call$data))) data <- as.data.frame(data)
  mf <- model.frame(formula, data)
  Terms <- attr(mf, "terms")
  
  X <- model.matrix(Terms, mf)
  y <- model.extract(mf, "response")
  
  ret <- ldm.default(y, X, ...)
  ret$call <- call
}
  
}

ldm_train <- function(y, X, gamma, lambda_1, lambda_2, cost)
{
  
  # construct kernel matrix
  m <- nrow(X)
  K <- matrix(nrow=m, ncol=m)
  for (i in 1:m) 
  {
    for (j in 1:m)
    {
      if (i==j) K[i,j] <- 1
      else if (i > j) K[i,j] <- K[j,i]
      else K[i,j] <- exp(-gamma * sum((X[i,] - X[j,])^2))
    }
  }
  
  beta <- rep(0, m)
  Q <- 4 * lambda_1 * ((m * t(K) %*% K) - (K %*% y %*% t(y) %*% t(K))) / (m^2) + K
  alpha <- lambda_2 / m * solve(Q) %*% K %*% y
  
  Y <- diag(y)
  A <- solve(Q) %*% K %*% Y
  h <- diag(Y %*% K %*% solve(Q) %*% K %*% Y)
  
  BETA_EPSILON <- 10^-6
  n <- 0
  while(TRUE) 
  {
    beta_old <- vector(length=m)
    del_beta <- vector(length=m)
    for(i in 1:m)
    {
      del_beta[i] <- (Y %*% K %*% alpha)[i] - 1
      beta_old[i] <- beta[i]
      beta[i] <- min(max(beta[i] - del_beta[i]/h[i], 0), cost)
      alpha <- alpha + (beta[i] - beta_old[i]) * A[,i]
    }
    if (sum((beta - beta_old)^2) < BETA_EPSILON) break
    if( 1000 < (n <- n + 1)) break
  }
  
  list(alpha = alpha, K=K)
}