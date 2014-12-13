predict.ldm <- function(object, newdata, type="response") {
  
  if(missing(newdata)) return(object$fitted)
  
  type <- c("response", "decision.values")[pmatch(type, c("response", "decision.values"))]
  
  ### ONLY WORKS FOR RBF KERNEL A.T.M. ###
  
  if(inherits(object, "ldm.formula")) {
    
    mf <- model.frame(object$formula, newdata)
    Terms <- attr(mf, "terms")
    attr(Terms, "intercept") <- 0
    
    newdata <- model.matrix(Terms, mf)
  }
  
  alpha <- object$alpha
  mod.matrix <- object$model.matrix
  d <- c(dim(newdata)[1], dim(mod.matrix))
  gamma <- object$gamma
  
  scaled <- object$scaled
  
  if (any(scaled)) {
    x.scale <- object$x.scale
    xtmp <- scale(newdata[ , scaled], center=x.scale[[1]], scale=x.scale[[2]])
    newdata[ , scaled] <- xtmp
  }
  rownames(newdata) <- NULL
  
  if (object$kernel != 0) {
    pred <- .C("ldmPredict", as.double(alpha), as.double(gamma), 
      as.double(newdata), as.integer(d), as.double(mod.matrix), 
      pred=double(d[1]), PACKAGE="rldm")$pred
  }
  else {
    pred <- newdata %*% t(mod.matrix) %*% as.matrix(alpha, ncol=1)
  }
  if (identical(type, "response")) pred <- as.factor(object$y.levels[sign(pred)/2 + 1.5])
  levels(pred) <- object$y.levels
  return(pred)
}
  
