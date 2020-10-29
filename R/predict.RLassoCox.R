predict.RLassoCox <- function(object, newx, ...){
    newx <- newx[,rownames(object$glmnetRes$beta)]  
    predict(object$glmnetRes, newx, ...)
}
