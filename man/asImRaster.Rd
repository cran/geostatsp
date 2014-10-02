\name{asImRaster}
\alias{asImRaster}


\title{Convert a raster to an im object}
\description{
Conversion between rasters and spatstat's im objects
}
\usage{
	asImRaster(X, ...)
}

\arguments{
	\item{X}{
A RasterLayer}
\item{...}{ additional arguments}
}

\details{
This function is not registered as an S3 method as doing so would require the spatstat package to be a dependency of geostatsp.
}

\examples{

myraster = raster(matrix(1:100,10,10),
	xmn=0,xmx=10,ymn=0,ymx=10)

if(require('spatstat',quietly=TRUE)){	
	myIm = asImRaster(myraster)
}
}


