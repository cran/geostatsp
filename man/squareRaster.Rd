\name{squareRaster-methods}
\docType{methods}
\alias{squareRaster}
\alias{squareRaster-methods}
\alias{squareRaster,matrix-method}
\alias{squareRaster,Raster-method}
\alias{squareRaster,Spatial-method}
\alias{squareRaster,Extent-method}



\title{Create a raster with square cells}
\description{
Given a raster object, an extent, or a bounding box, a raster of with square cells and having the 
extent and number of cells specified is returned.
}
\usage{
\S4method{squareRaster}{matrix}(x,cells=NULL, buffer=0)
\S4method{squareRaster}{Raster}(x,cells=NULL, buffer=0)
\S4method{squareRaster}{Spatial}(x,cells=NULL, buffer=0)
\S4method{squareRaster}{Extent}(x,cells=NULL, buffer=0)
}

\arguments{
	\item{x}{
	A bounding box from a  
	\code{SpatialPoints} or \code{SpatialPolygons} object or an 
	\code{Extent} from a Raster. }
\item{cells}{The number of cells in the x direction.  If NULL the number of columns of x is used.  }
\item{buffer}{Additional area to add around the resulting raster}
}

\value{
A \code{rasterLayer} with square cells
}	

\examples{

myraster = raster(matrix(0,10,10),xmn=0,xmx=10,ymn=0,ymx=12.3)

squareRaster(myraster)

squareRaster(extent(myraster), cells=10)

squareRaster(bbox(myraster), cells=10)

}