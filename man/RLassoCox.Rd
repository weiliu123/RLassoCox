\name{RLassoCox}
\alias{RLassoCox}
\title{
Reweighted Lasso-Cox model
}
\description{
A rewighted Lasso-Cox model for survival prediction and biomarker discovery.
}
\usage{
RLassoCox(x, y, globalGraph = NULL, Gamma = 0.3, DEBUG = TRUE, 
            standardize = TRUE, ...)
}
\arguments{
    \item{x}{
a n x p matrix of gene expression measurements with n samples and p genes.
}
    \item{y}{
a n x 2 matrix of survival data. The two columns represent disease status 
'status' and survival time 'time' respectively. 
}
    \item{globalGraph}{
An \code{igraph} R object containing the interaction network.
}
    \item{Gamma}{
A numeric value. The restart probability in directed random walk. Default 
is \code{Gamma = 0.3}.
}
    \item{DEBUG}{
Logical. Should debugging information be shown.
}
    \item{standardize}{
Logical flag for \code{x} standardization, prior to fitting the model. 
Default is \code{TRUE}.
}
    \item{\dots}{
Arguments to be passed to \code{glmnet} in R package \code{glmnet}.
}
}
\details{
RLassoCox integrates gene interaction information into the Lasso-Cox model 
for accurate survial prediction and biomarker discovery.
}
\value{
    \item{glmnetRes}{An object of class "glmnet"}
    \item{PT}{The topological weights of genes}
}
\references{
Integration of gene interaction information into a reweighted Lasso-Cox model 
for accurate survival prediction. To be published.
}
\author{
Wei Liu
}


\seealso{
    \code{\link{predict}}
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
}
