\name{predict.cvRLassoCox}
\alias{predict.cvRLassoCox}
\title{
Make predictions from a cross-validated RLasso-Cox model
}
\description{
This function makes predictions from a cross-validated RLasso-Cox model, using 
the optimal value chosen for lambda.
}
\usage{
\method{predict}{cvRLassoCox}(object, newx, ...)
}
\arguments{
    \item{object}{
cross-validated RLasso-Cox model
}
    \item{newx}{
A matrix with new samples to predict.
}
    \item{\dots}{
Arguments to be passed to \code{predict.cv.glmnet} in R package \code{glmnet}.
}
}
\value{
Predicted results of new patients in \code{newx}.
}
\examples{
library("survival")
library("igraph")
library("glmnet")
library("Matrix")

data(dGMMirGraph)
data(mRNA_matrix)
data(survData)

trainSmpl.Idx <- sample(1:dim(mRNA_matrix)[1], floor(2/3*dim(mRNA_matrix)[1]))
testSmpl.Idx <- setdiff(1:dim(mRNA_matrix)[1], trainSmpl.Idx)
trainSmpl <- mRNA_matrix[trainSmpl.Idx ,]
testSmpl <- mRNA_matrix[testSmpl.Idx ,]

cv.res <- cvRLassoCox(x=trainSmpl, y=survData[trainSmpl.Idx ,], 
                        globalGraph=dGMMirGraph, nfolds = 5)
lp <- predict(object = cv.res, newx = testSmpl, 
                            s = "lambda.min")
}
