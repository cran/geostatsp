\name{openmap}
\alias{openmap}
\alias{openmap.default}
\alias{openmap.matrix}
\alias{openmap.Extent}

\title{Get map from a server}
\description{
	This function is a wrapper for the 
	\code{\link[OpenStreetMap]{openmap}} 
function in the \link[OpenStreetMap]{openmap} package, where
the limits can be passed as a bounding box or an Extent.
}
\usage{
	openmap(upperLeft,  ...)
	\method{openmap}{matrix}(upperLeft,...)
	\method{openmap}{Extent}(upperLeft,...)
	}

\arguments{
	\item{upperLeft}{
	A bounding box from a  
	SpatialPoints or SpatialPolygons object or an 
	Extent from a Raster.  The name 'upperLeft' is retained for compatibility with 
	\code{\link[OpenStreetMap]{openmap}}, but this argument specifies a bounding box and not 
	an upper left corner. }
\item{...}{ additional arguments, 
see \code{\link[OpenStreetMap]{openmap}} in the OpenStreetMap package.}
}

\value{
	an object of class \code{OpenStreetMap}, in a lat-long projection.
}	

\examples{

myraster = raster(matrix(0,10,10),xmn=0,xmx=10,ymn=0,ymx=10)

plot(myraster)

# only do the following if the the OpenStreetMap.org web site can be accessed
#if(length(nsl("www.OpenStreetMap.org"))) {
	mytiles = openmap(myraster@extent, zoom=4,  type="osm")

	plot(myraster)
	plot(mytiles, add=TRUE)

#	mypoints = SpatialPoints(10*cbind(runif(10), runif(20)))

#	mytiles=openmap(mypoints@bbox, zoom=4,  type="bing")

#	plot(mypoints)
#	plot(mytiles,add=T)
#	points(mypoints, col='yellow')
#}

}