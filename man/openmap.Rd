\name{openmap}
\alias{openmap}
\alias{getTiles}
\title{Get map from a server}
\description{
	This function is a wrapper for the 
	\code{\link[OpenStreetMap]{openmap}} 
function in the \link[OpenStreetMap]{openmap} package
}
\usage{
	openmap(x, zoom, 
	path="http://tile.openstreetmap.org/",
	maxTiles = 9,
	crs, cacheDir=tempdir(), verbose=FALSE)
	getTiles(xlim,ylim,zoom,path,maxTiles = 16,cacheDir=tempdir(),
		timeOut=5*24*60*60,verbose=FALSE)
}

\arguments{
	\item{x}{
	Either	 raster, SPDF, matrix corresponding to a bounding box or an Extent.  }
\item{zoom	}{
the zoom level, when missing it will be determined by maxTiles.}
\item{path}{	
the tile server from which to get the map, see \url{http://wiki.openstreetmap.org/wiki/Tiles#Servers}.}
\item{
maxTiles}{	
If zoom is null, zoom will be chosen such that the 
number of map tiles is less than or equl to this number.}
\item{crs}{Projection for the output, defaulting to the same projection as
\code{x}.  If \code{x} has no projection, 
for instance when \code{x} is a matrix or extent, 
\code{crs} is also used as the projection of \code{x}.  If \code{crs} is missing
and \code{x} has no crs, long-lat is used.}

  \item{cacheDir}{
    A directory for caching tiles retrieved over the web. Not used if
    tiles are coming from a local directory.
  }
\item{xlim}{
    Length-2 vector of min and max X coordinate (latitude)
  }
  \item{ylim}{
    Length-2 vector of min and max Y coordinate (longitude)
  }

  \item{timeOut}{
    Tiles in the cache directory older than this (in minutes) are
    re-retrieved.
  }
  \item{verbose}{
    If TRUE, give information on where tiles are coming from, cache
    hits, etc
  }
}

\details{
These functions download, display, and manipulate map tiles stored in a
standard way either on a web server or a local folder.

Map tiles are a set of PNG images that span the world at a set of zoom
levels. Zoom level 0 has one tile that covers the whole world at a
resolution of 256x256 pixels. Zoom level 1 has four 256x256 pixel tiles
in a 2x2 pattern over the whole world. In general, zoom level n has 2^n
by 2^n tiles. Zoom levels go up to about 17 or 18 depending on the tile
server.

The function getTiles() retrieves tiles to cover a given latitude and
longitude range at a given zoom level. Note that since tiles are
generated at a fixed set of lat-long ranges the returned tiles will
overlap the requested area. There may be some way of clipping to the
requested area in a future revision.

The getTilePaths() function returns information about tiles in a region
at a zoom level. The tile source path is given but tiles are not
requested. A list is returned with one element per tile giving the path,
the I and J coordinates, the zoom level, and the geographic bounds for
the tile. 

The tileBbox() function returns the bounds of the tiles by simply
constructing the bounding box of all the constitutent
SpatialGridDataFrame objects.

http://wiki.openstreetmap.org/index.php/Tile_usage_policy

}


\value{
	\code{openmap} returns a RasterBrick 
	\code{\link[raster]{brick}}, with 'red', 'green' and 'blue' layers. \code{getTiles} returns a list
	of RasterBricks
}	



\examples{


myraster = raster(matrix(0,10,10),xmn=8,xmx=18,ymn=0,ymx=10, crs="+proj=longlat")
values(myraster) = seq(0,1,len=ncell(myraster))

myPoints = SpatialPoints(myraster, proj4string=CRS(proj4string(myraster)))[
	seq(1,ncell(myraster),len=5)]




\dontrun{
	par(mar=c(0,0,0,0))
	mytiles = openmap(myraster)

	plot(myraster)
	plot(mytiles, add=TRUE)
	points(myPoints,col='red')



	mytiles = openmap(myPoints,  
		path="b.www.toolserver.org/tiles/bw-mapnik", verbose=TRUE)
	plot(myPoints)
	plot(mytiles, add=TRUE)
	points(myPoints,col='red')
	
}

}