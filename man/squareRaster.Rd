\name{squareRaster-methods}
\docType{methods}
\alias{squareRaster}
\alias{squareRaster-methods}
\alias{squareRaster,matrix-method}
\alias{squareRaster,SpatRaster-method}
\alias{squareRaster,SpatVector-method}
\alias{squareRaster,SpatExtent-method}



\title{Create a raster with square cells}
\description{
Given a raster object, an extent, or a bounding box, a raster of with square cells and having the 
extent and number of cells specified is returned.
}
\usage{
\S4method{squareRaster}{matrix}(x,cells=NULL, buffer=0)
\S4method{squareRaster}{SpatRaster}(x,cells=NULL, buffer=0)
\S4method{squareRaster}{SpatVector}(x,cells=NULL, buffer=0)
\S4method{squareRaster}{SpatExtent}(x,cells=NULL, buffer=0)
}

\arguments{
	\item{x}{
	A spatial object}
\item{cells}{The number of cells in the x direction.  If NULL the number of columns of x is used.  }
\item{buffer}{Additional area to add around the resulting raster}
}

\value{
A \code{SpatRaster} with square cells
}	

\examples{

myraster = rast(matrix(0,10,10),extent=c(0,10,0,12.3))

squareRaster(myraster)

squareRaster(myraster, buffer=3, cells=20)

squareRaster(ext(myraster), cells=10)


}