\name{as.im}
\alias{as.im}
\alias{as.im.RasterLayer}
\alias{as.im.default}


\title{Convert a raster to an im object}
\description{
Conversion between rasters and spatstat's im objects
}
\usage{
	as.im(X, ...)
	\method{as.im}{RasterLayer}(X,...)
	\method{as.im}{default}(X,...)
	}

\arguments{
	\item{X}{
A RasterLayer or an object acceptable to 
\code{\link[spatstat]{as.im}}.}
\item{...}{ additional arguments}
}



\examples{

myraster = raster(matrix(1:100,10,10),
	xmn=0,xmx=10,ymn=0,ymx=10)
\dontrun{
myIm = as.im(myraster)
}
}


