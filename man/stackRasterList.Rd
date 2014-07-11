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
stackRasterList(x, template = x[[1]], method = "ngb", mc.cores=NULL)
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
\item{mc.cores}{If non-null, \code{\link[parallel]{mcmapply}} is used with this argument specifying the number of cores.}
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



mylist = list(a=raster(matrix(1:36, 6, 6,byrow=TRUE), 0,1000,0,1000, 
       crs="+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs"),
	b=raster(matrix(1:144, 12, 12), -85.49, -85.48, 0, 0.01,
	   crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "),
	c=raster(matrix(1:100, 10, 10), -5000,5000,-5000,5000, 
       crs="+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
	)
	
if(.Platform$OS.type=='unix') {	
	mystack = stackRasterList(mylist,mc.cores=2)
	mystack
}

plot(mystack[["b"]], main="stack b")
plot(mystack[['a']],add=TRUE,col=grey(seq(0,1,len=12)),alpha=0.8,legend=FALSE)

}

