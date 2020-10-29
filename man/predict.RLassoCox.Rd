\name{predict.RLassoCox}
\alias{predict.RLassoCox}
\title{
Make predictions from a RLasso-Cox model
}
\description{
This function predicts the risk of new samples from a fitted \code{RLasso-Cox}
model.
}
\usage{
\method{predict}{RLassoCox}(object, newx, ...)
}
\arguments{
    \item{object}{
Fitted \code{"RLassoCox"} model object.
}
    \item{newx}{
A matrix with new samples to predict.
}
    \item{\dots}{
Arguments to be passed to \code{predict.glmnet} in R package \code{glmnet}.
}
}
\value{
Predicted results of new patients in \code{newx}.
}
\author{
Wei Liu
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

res <- RLassoCox(x=trainSmpl, y=survData[trainSmpl.Idx ,],
                globalGraph=dGMMirGraph)
lp <- predict(object = res, newx = testSmpl)
}

