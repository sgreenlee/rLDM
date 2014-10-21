
ldm <- function(x, ...) UseMethod("ldm")

ldm.default <- function(y, X, kernel = 'radial', cost = 10, gamma=1, degree=2, lambda_1=1/4, lambda_2=1/4)
{
  kernel_type <- pmatch(kernel, c('linear', 'polynomial', 'radial', 'sigmoid')) - 1
  
  if (!is.factor(y)) stop("y must be a factor")
  if (length(levels(y)) != 2) stop("y must be a two-level factor: multi-class not implemented yet")
  else y <- c(-1,1)[unclass(y)]
  
  X <- scale(X)
  
  tuning_params <- c(cost, lambda_1, lambda_2, degree, gamma)
  n <- dim(X)[1]
  m <- dim(X)[2]
  
  gmat <- .C("makeGramMatrix", as.double(X), as.integer(n), as.integer(m),
              G=double(n^2), as.integer(kernel_type), as.double(tuning_params),
             PACKAGE="rldm")
  G <- matrix(gmat$G, nrow=n)
  Gy <- G %*% y
  Y <- diag(y)
  GY <- G %*% Y
  Q <- 4 * lambda_2 * ((t(G) %*% G) / n - Gy %*% t(Gy) / n^2) + G
  invQ <- solve(Q)
  A <- invQ %*% GY
  h <- diag(t(GY) %*% A)
  
  alpha <- lambda_2 / n * invQ %*% Gy 
  beta <- double(n)
  beta_old <- double(n)
  del_beta <- double(n)
  
  cret <- .C("coordDescent", as.integer(n), as.double(A), as.double(GY), 
             alpha = as.double(alpha), beta = as.double(beta), as.double(beta_old), 
             as.double(del_beta), as.double(h), as.double(cost), PACKAGE="rldm")
  
  ret <- list(alpha=cret$alpha)
  ret$G <- G
  ret$fitted <- sign(G %*% ret$alpha)
  ret$cost <- cost
  ret$lambda_1 <- lambda_1
  ret$lambda_2 <- lambda_2
  ret$gamma <- gamma
  ret$degree <- degree
  class(ret) <- 'ldm'
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