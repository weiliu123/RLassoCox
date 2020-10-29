\name{mRNA_matrix}
\alias{mRNA_matrix}
\docType{data}
\title{
The expression data
}
\description{
An example of GBM expression data. We acknowledge the TCGA Research Network 
for generating the GBM datasets.
}
\usage{data("mRNA_matrix")}
\format{
    The format is:
    num [1:314, 1:4853] 0.562167 0.022435 -0.000102 -0.719444 0.620269 ...
    - attr(*, "dimnames")=List of 2
    ..$ : chr [1:314] "TCGA-02-0001" "TCGA-02-0003" "TCGA-02-0006" 
    ..$ : chr [1:4853] "90993" "4313" "26248" "57680" ...
}
\examples{
data(mRNA_matrix)
}
\keyword{datasets}
