predict.ldm <- function(object, newdata, type=c("response", "decision.values")) {
  
  if(missing(newdata)) return(object$fitted)
  
  ### ONLY WORKS FOR RBF KERNEL A.T.M. ###
  
  if(inherits(object, "ldm.formula")) {
    
    # construct model matrix from formula
  
  }
  
  ####                ####
  #       SCALE DATA     #
  ####                ####
  
  alpha <- object$alpha
  mod.matrix <- object$model.matrix
  d <- c(dim(newdata)[1], dim(mod.matrix))
  gamma <- object$gamma
  
  pred <- .C("ldmPredict", as.double(alpha), as.double(gamma), 
             as.double(newdata), as.integer(d), as.double(mod.matrix), 
             pred=double(d[1]), PACKAGE="rldm")$pred
  if (identical(type, "response")) pred <- as.factor(sign(pred))
  return(pred)
}
  
