\name{pcPriorRange}
\alias{pcPriorRange}
\alias{pcPrior}
\title{
PC prior for range parameter
}
\description{
Creates a penalized complexity prior for the range parameter
}
\usage{
pcPriorRange(q, p=0.5, cellSize=1)
}
 
\arguments{
\item{q}{Lower quantile for the range parameter}
\item{p}{probability that the range is below this quantile, defaults to the median}
\item{cellSize}{size of grid cells, can be a raster.}
}
\details{
q is the quantile in spatial units, usually meters, and the scale parameter follows an exponential distribution. 	A prior PC prior distribution for the range parameter in units of grid cells, which INLA requires, is computed.  
}
\value{
	A list with
	\item{lambda}{parameter for the exponential distribution (for scale in units of cells), in the same parametrization as dexp}
	\item{priorScale}{matrix with x and y columns with prior of scale parameter}
	\item{priorRange}{matris with x and y columns with prior of range parameter, in meters (or original spatial units)}
	\item{inla}{character string specifying this prior in inla's format}
}


\examples{
# pr(range < 100km) = 0.1, 200m grid cells 
	x = pcPriorRange(q=100*1000, p=0.1, cellSize = 200)
	plot(x$priorRange \%*\% diag(c(1/1000000, 1000000)), 
	  type='l', xlab="range, 1000's km", ylab='dens',
 	  xlim = c(0, 1500/1000))
	cat(x$inla)
}