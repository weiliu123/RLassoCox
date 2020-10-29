predict.cvRLassoCox <- function(object, newx,...){
    newx <- newx[,rownames(object$glmnetRes$glmnet.fit$beta)]
    predict(object$glmnetRes, newx, ...)
}
