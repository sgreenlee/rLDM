predict.ldm <- function(object, newdata, type="response") {
  
  if(missing(newdata)) return(object$fitted)
  
  type <- c("response", "decision.values")[pmatch(type, c("response", "decision.values"))]
  
  ### ONLY WORKS FOR RBF KERNEL A.T.M. ###
  
  if(inherits(object, "ldm.formula")) {
    
    # construct model matrix from formula
  
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
  
  pred <- .C("ldmPredict", as.double(alpha), as.double(gamma), 
             as.double(newdata), as.integer(d), as.double(mod.matrix), 
             pred=double(d[1]), PACKAGE="rldm")$pred
  if (identical(type, "response")) pred <- as.factor(sign(pred))
  return(pred)
}
  
