\name{rw}
\alias{rw}
\title{
Directed Random Walk
}
\description{
The directed random walk algorithm proposed by Liu et al (2013).
}
\usage{
rw(W, p0, gamma)
}
\arguments{
    \item{W}{
The adjacency matrix of the gene interaction network.
}
    \item{p0}{
A vector containing the initial weights of genes in the gene interaction 
network.
}
    \item{gamma}{
A numeric value. The restart probability in directed random walk.
}
}
\details{
This function implements the directed random walk algorithm proposed by 
Liu et al (2013). It evaluates the topological weight of each gene according 
to its topological importance in the gene interaction network. The genes 
that close to many other genes that have large weights will receive larger 
weights. The final weights reflect the topological importances of genes in 
the gene interaction network. 
}
\value{
A matrix containing the topological weights of nodes in \code{igraphM}.
}
\references{
Liu, W., et al., Topologically inferring risk-active pathways toward precise 
cancer classification by directed random walk. Bioinformatics, 2013. 29(17): 
p. 2169-77.
}
\author{
Wei Liu <freelw@qq.com>
}
