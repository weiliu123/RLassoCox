\name{RLassoCox-package}
\alias{RLassoCox-package}
\docType{package}
\title{
\packageTitle{RLassoCox}
}
\description{
\packageDescription{RLassoCox}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{RLassoCox}
\packageIndices{RLassoCox}
Very simple to use. Accepts x,y data for the RLasso-Cox model, and makes 
predictions for new samples. 

\code{RLassoCox} A rewighted Lasso-Cox model for survival prediction and 
biomarker discovery.
\code{predict.RLassoCox} This function predicts the risk of new samples 
from a fitted \code{RLasso-Cox} model.
\code{cvRLassoCox} Does k-fold cross-validation for the RLasso-Cox model, 
produces a plot, and returns a value for lambda
\code{predict.cvRLassoCox} This function makes predictions from a 
cross-validated RLasso-Cox model, using the optimal value chosen for lambda.
}
\author{
\packageAuthor{RLassoCox}

Maintainer: \packageMaintainer{RLassoCox}
}
\references{
Integration of gene interaction information into a reweighted Lasso-Cox model 
for accurate survival prediction. To be published.
}
\keyword{ package }
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

res <- RLassoCox(x=trainSmpl, y=survData[trainSmpl.Idx ,], 
                globalGraph=dGMMirGraph)
lp <- predict(object = res, newx = testSmpl)

cv.res <- cvRLassoCox(x=trainSmpl, y=survData[trainSmpl.Idx ,], 
                        globalGraph=dGMMirGraph, nfolds = 5)
cv.lp <- predict(object = cv.res, newx = testSmpl, 
                            s = "lambda.min")
}
