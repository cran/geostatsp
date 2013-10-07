\name{squareRaster}
\alias{squareRaster}
\alias{squareRaster.BasicRaster}
\alias{squareRaster.RasterLayer}
\alias{squareRaster.matrix}
\alias{squareRaster.Extent}
\alias{squareRaster.SpatialPoints}
\alias{squareRaster.SpatialPointsDataFrame}
\alias{squareRaster.SpatialPolygons}
\alias{squareRaster.SpatialPolygonsDataFrame}

\title{Create a raster with square cells}
\description{
Given a raster object, an extent, or a bounding box, a raster of with square cells and having the 
extent and number of cells specified is returned.
}
\usage{
	squareRaster(x,  cells=NULL)
	\method{squareRaster}{RasterLayer}(x,cells=NULL)
	\method{squareRaster}{BasicRaster}(x,cells=NULL)
	\method{squareRaster}{matrix}(x,cells=NULL)
	\method{squareRaster}{Extent}(x,cells=NULL)
	\method{squareRaster}{SpatialPoints}(x,cells=NULL)
	\method{squareRaster}{SpatialPointsDataFrame}(x,cells=NULL)
	\method{squareRaster}{SpatialPolygons}(x,cells=NULL)
	\method{squareRaster}{SpatialPolygonsDataFrame}(x,cells=NULL)
	}

\arguments{
	\item{x}{
	A bounding box from a  
	SpatialPoints or SpatialPolygons object or an 
	Extent from a Raster. }
\item{cells}{The number of cells in the x direction.  If NULL the number of columns of x is used.  }
}

\value{
A raster.
}	

\examples{

myraster = raster(matrix(0,10,10),xmn=0,xmx=10,ymn=0,ymx=12.3)

squareRaster(myraster)

squareRaster(extent(myraster), cells=10)

squareRaster(bbox(myraster), cells=10)

}