\name{nn32}
\alias{nn32}
\alias{nn64}
\alias{nn128}
 
\docType{data}
\title{
Nearest neighbour matrices
}
\description{
Some large nearest neighbour matrices.
}
\usage{data(nn32)}
\format{
Sparse matrices produced by \code{\link{NNmat}}, showing the neighbours of cells on a 64 by 64 grid 
in \code{nn64},  a 128 by 128 grid in \code{nn128}, and 256 in \code{nn256}.
}
 

\examples{
data(nn32)
dim(nn32)
dim(nn64)
dim(nn128)
x = 50
y=31
matrix(nn128[,128*(y-1)+x], 128, 128)[ x+seq(-4,4),y+seq(-4, 4)]

length(nn128@x)
}


\keyword{datasets}
