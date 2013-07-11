\name{stackRasterList}
\alias{stackRasterList}

\title{
Converts a list of rasters, possibly with different projections and resolutions, to a single raster stack.
}
\description{
This function is intended to be used prior to passing covariates to \link{krige} in order
for the rasters for all covariates to have the same projection and same resolution.  

}
\usage{
stackRasterList(x, template = x[[1]], method = "ngb")
}
\arguments{
  \item{x}{
A list of rasters
}
  \item{template}{
A raster whose projection and resolution all other rasters will be aligned with.  Defaults to the first raster in \code{x}
}
  \item{method}{
The method to use, either "ngb", or "bilinear".  Can be a vector of the same length as x to specify different methods for each raster.  If \code{method} has names 
which correspond to the names of \code{x}, the names will be used instead of the order to assign methods to rasters.
}
}

\value{
A raster stack.
}


\examples{

mylist = list(a=raster(matrix(1:9, 3, 3), 0,1,0,1, 
       crs="+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs"),
	b=raster(matrix(1:25, 5, 5), -1, 2, -1, 2,
	   crs="+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
	)
mystack = stackRasterList(mylist)
mystack



mylist = list(a=raster(matrix(1:9, 3, 3), 0,1000,0,1000, 
       crs="+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs"),
	b=raster(matrix(1:25, 5, 5), -85.49, -85.48, 0.003, 0.009,
	   crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
	)
	
mystack = stackRasterList(mylist)
mystack


}

