\name{dGMMirGraph}
\alias{dGMMirGraph}
\docType{data}
\title{
The KEGG network
}
\description{
The KEGG network constructed by the R package \code{iSubpathwayMiner}. 
}
\usage{data("dGMMirGraph")}
\format{
An \code{igraph} R object.
}
\details{
There are 7159 nodes and 39930 edges in \code{dGMMirGraph}. Each node in 
the graph represents a gene/miRNA/metabolite. The KEGG network is used to 
evaluate the topological importance of genes by the random walk method.
}
\examples{
data(dGMMirGraph)
}
\keyword{datasets}
