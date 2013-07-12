\name{excProb}
\alias{excProb}
 

\title{
 Exceedance probabilities
}
\description{
Calculate exceedance probabilities pr(X > threshold) from an INLA result.}
\usage{
excProb(marginals, threshold=0, template=NULL, templateIdCol=NULL)
}
 
\arguments{
\item{marginals}{
A list of two-column matrices with columns named \code{x} and \code{y} containing the 
posterior distributions of random effects, as produced by \code{inla}.
}
\item{threshold}{ the value which 
}
\item{template}{A \code{Raster} or \code{SpatialPolygonsDataFrame} or \code{SpatialPointsDataFrame} 
object which the results will be contained in.
}
\item{templateIdCol}{The data column in \code{template} corresponding to names of \code{marginals}
}
}
\details{
	Numerically integrates the marginal posteriors to obtain pr(X > threshold).
}
\value{
Either a vector of exceedance probabilities or an object of the same class as \code{template}.
}

\examples{
\dontrun{
	swissFit =  glgm(swissRain, cells=40, formula="lograin",
	covariates=swissAltitude)
	swissExc = excProb(swissFit$inla$marginals.random$space, log(300),
		swissRain)
}

}


 	