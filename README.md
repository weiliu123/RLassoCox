# RLassoCox
A reweighted Lasso-Cox model for survival prediction and biomarker discovery.

# Details
Package: RLassoCox

Type: Package

Title: A reweighted Lasso-Cox by integrating gene interaction information

Version: 0.99.3

Date: 2020-11-20

Authors@R: c(person(given = "Wei", family = "Liu", email = "freelw@qq.com", role = c("cre", "aut"),comment = c(ORCID = "0000-0002-5496-3641")))

Depends: R (>= 4.1), glmnet

Imports: Matrix, igraph, survival, stats

Description: RLassoCox is a package that implements the RLasso-Cox model proposed by Wei Liu. The RLasso-Cox model integrates gene interaction information into the Lasso-Cox model for accurate survival prediction and survival biomarker discovery. It is based on the hypothesis that topologically important genes in the gene interaction network tend to have stable expression changes. The RLasso-Cox model uses random walk to evaluate the topological weight of genes, and then highlights topologically important genes to improve the generalization ability of the Lasso-Cox model. The RLasso-Cox model has the advantage of identifying small gene sets with high prognostic performance on independent datasets, which may play an important role in identifying robust survival biomarkers for various cancer types.

License: Artistic-2.0

biocViews: Survival, Regression, GeneExpression, GenePrediction, Network

BugReports: https://github.com/weiliu123/RLassoCox

BiocType: Software

Suggests: knitr

VignetteBuilder: knitr

# Index of help topics:
RLassoCox Reweighted Lasso-Cox model

cvRLassoCox Cross-validation for the RLasso-Cox model

dGMMirGraph The KEGG network

g.HuRI.EntrezID The HuRI network

mRNA_matrix The expression data

predict.RLassoCox Make predictions from a RLasso-Cox model

predict.cvRLassoCox Make predictions from a cross-validated RLasso-Cox model

rw Directed Random Walk

survData Survival data

# Examples
data(dGMMirGraph)

data(mRNA_matrix)

data(survData)

trainSmpl.Idx <- sample(1:dim(mRNA_matrix)[1], floor(2/3*dim(mRNA_matrix)[1]))

testSmpl.Idx <- setdiff(1:dim(mRNA_matrix)[1], trainSmpl.Idx)

trainSmpl <- mRNA_matrix[trainSmpl.Idx ,]

testSmpl <- mRNA_matrix[testSmpl.Idx ,]

res <- RLassoCox(x=trainSmpl, y=survData[trainSmpl.Idx ,], globalGraph=dGMMirGraph) 

lp <- predict(object = res, newx = testSmpl)

cv.res <- cvRLassoCox(x=trainSmpl, y=survData[trainSmpl.Idx ,], globalGraph=dGMMirGraph, nfolds = 5) 

cv.lp <- predict(object = cv.res, newx = testSmpl, s = "lambda.min")
