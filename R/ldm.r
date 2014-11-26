
ldm <- function(x, ...) UseMethod("ldm")

ldm.default <- function(y, X, kernel = 'radial', cost = 10, gamma=1, degree=2, 
                        a=1, c=1, r=0, lambda_1=1/4, lambda_2=1/4, scale=T)
{
  kernel_type <- pmatch(kernel, c('linear', 'polynomial', 'radial', 'sigmoid')) - 1
  
  if (!is.factor(y)) stop("y must be a factor")
  if (length(levels(y)) != 2) stop("y must be a two-level factor: multi-class not implemented yet")
  y.levels <- levels(y)
  y <- c(-1,1)[unclass(y)]
  
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
  
  tuning_params <- c(r, a, c, degree, gamma)
  n <- dim(X)[1]
  m <- dim(X)[2]
  
  if (kernel_type != 0) {
    gmat <- .C("makeGramMatrix", as.double(X), as.integer(n), as.integer(m),
              G=double(n^2), as.integer(kernel_type), as.double(tuning_params),
             PACKAGE="rldm")
  }
  G <- matrix(gmat$G, ncol=n)
  GY <- matrix(gmat$G, nrow=n) %*% diag(y)
  Q <- 4 * lambda_1 * ((t(G) %*% G) / n - GY %*% matrix(1, ncol=n, nrow= n) %*% t(GY) / n^2) + G
  A <- solve(Q, GY)
  alpha <- lambda_2 / n * A %*% matrix(rep(1, times=n), ncol=1)
  h <- diag(t(GY) %*% A)
  
  beta <- double(n)
  beta_old <- double(n)
  del_beta <- double(n)
  
  cret <- .C("coordDescent", as.integer(n), as.double(A), as.double(GY), 
             alpha = as.double(alpha), beta = as.double(beta), as.double(beta_old), 
             as.double(del_beta), as.double(h), as.double(cost), PACKAGE="rldm")
  
  ret <- list(alpha=cret$alpha)
  ret$G <- G
  ret$fitted <- y.levels[sign(G %*% ret$alpha)/2 + 1.5]
  ret$y.levels <- y.levels
  ret$cost <- cost
  ret$lambda_1 <- lambda_1
  ret$lambda_2 <- lambda_2
  ret$gamma <- gamma
  ret$degree <- degree
  ret$model.matrix <- X
  ret$scaled <- scale
  ret$x.scale <- x.scale
  class(ret) <- 'ldm'
  return(ret)
}
  
ldm.formula <- function(formula, data = NULL, ... )
{
  call <- match.call()
  if(is.matrix(eval.parent(call$data))) data <- as.data.frame(data)
  mf <- model.frame(formula, data)
  Terms <- attr(mf, "terms")
  attr(Terms, "intercept") <- 0
  
  X <- model.matrix(Terms, mf)
  y <- model.extract(mf, "response")
  
  ret <- ldm.default(y, X, ...)
  ret$call <- call
  ret$formula <- formula
  class(ret) <- c('ldm.formula', class(ret))
  return(ret)
}