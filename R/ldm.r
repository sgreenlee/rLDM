
ldm <- function(x, ...) UseMethod("ldm")

ldm.default <- function(y, X, kernel = 'radial', method = 'cd', cost = 10, gamma=1,
                        degree=2, a=1, c=1, r=0, lambda_1=1/4, lambda_2=1/4, eta = 0.01,
                        scale=T)
{
  library(MASS)
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
  
  method <- c('cd', 'asgd')[pmatch(method, c('cd', 'asgd'))]
  if(identical(method, 'asgd')){
    w <- double(m)
    w_bar <- double(m)
    gradient <- double(m)
    cret <- .C("avSGD", as.integer(n), as.integer(m), as.double(X), as.double(y), 
               as.double(cost), as.double(eta), as.double(w), as.double(w_bar), 
               as.double(gradient), as.double(lambda_1), as.double(lambda_2), 
               package='rldm')
    ret <- list(alpha=NULL)
    ret$G <- NULL
    ret$method <- 'ASGD'
    ret$fitted <- y.levels[sign(G %*% ret$alpha)/2 + 1.5]
    ret$y.levels <- y.levels
    ret$cost <- cost
    ret$eta <- eta
    ret$w <- cret$w_bar
    ret$lambda_1 <- lambda_1
    ret$lambda_2 <- lambda_2
    ret$gamma <- NULL
    ret$degree <- NULL
    ret$model.matrix <- X
    ret$scaled <- scale
    ret$x.scale <- x.scale
    class(ret) <- 'ldm'
    return(ret)
  }
  
  if (kernel_type != 0) {
    gmat <- .C("makeGramMatrix", as.double(X), as.integer(n), as.integer(m),
              G=double(n^2), as.integer(kernel_type), as.double(tuning_params),
             PACKAGE="rldm")
    G <- matrix(gmat$G, ncol=n)
  }
  else {
    G <- X %*% t(X)
  }
  GY <- G %*% diag(y)
  Q <- 4 * lambda_1 * ((t(G) %*% G) / n - GY %*% matrix(1, ncol=n, nrow= n) %*% t(GY) / n^2) + G
  Qinv <- ginv(Q)
  A <- Qinv %*% GY
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
  ret$method <- 'CD'
  ret$fitted <- y.levels[sign(G %*% ret$alpha)/2 + 1.5]
  ret$y.levels <- y.levels
  ret$cost <- cost
  ret$eta <- NULL
  ret$w <- NULL
  ret$kernel <- kernel_type
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