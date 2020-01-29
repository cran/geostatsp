\name{conditionalGmrf}
\alias{conditionalGmrf}
 \title{
Conditional distribution of GMRF
}
\description{
Distribution of Gaussian Markov Random Field conditional on data observed with noise on the same grid.
}
\usage{
conditionalGmrf(param, Yvec, Xmat, NN, 
	template = NULL, mc.cores = 1, 
	cellsPerLoop = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
vector of named parameters
}
  \item{Yvec}{
vector of observed data, or matrix with each column being a realisation.
}
  \item{Xmat}{
Matrix of covariates.
}
  \item{NN}{
 nearest neighbour matrix
}
  \item{template}{
Raster on which the GMRF is defined
}
  \item{mc.cores}{
passed to \code{\link[parallel:mclapply]{mcmapply}}
}
  \item{cellsPerLoop}{
number of cells to compute simultaneously.  Larger values consume more memory but 
result in faster computation.
}
  \item{\dots}{
  additional arguments passed to \code{\link{maternGmrfPrec}}
}
}
\value{
Raster image with layers containing conditional mean and standard deviation.
}
\author{
Patrick Brown
}

\seealso{
\code{\link{maternGmrfPrec}}, \code{\link{lgm}}
}

\keyword{models}

