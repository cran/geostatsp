\name{GNcities}
\alias{GNcities}
\alias{GNcities.default}
\alias{GNcities.matrix}
\alias{GNcities.Extent}

\title{Retrieve city names and locations}
\description{
	This function is a wrapper for the 
	\code{\link[geonames]{GNcities}} 
function in the \link[geonames]{geonames} package, where
the limits can be passed as a bounding box or an Extent.
}
\usage{
	GNcities(north,  ...)
	\method{GNcities}{matrix}(north,...)
	\method{GNcities}{Extent}(north,...)
	}

\arguments{
	\item{north}{
	A bounding box from a  
	SpatialPoints or SpatialPolygons object or an 
	Extent from a Raster.}
\item{...}{ additional arguments, 
see \code{\link[geonames]{GNcities}} in the geonames package.}
}

\value{A SpatialPointsDataFrame}


\examples{

mybbox = rbind(c(5,10),c(45,50))

\dontrun{

cities=GNcities(mybbox, max=5)

plot(cities)
text(cities, labels=cities$name, col='red')

}


}