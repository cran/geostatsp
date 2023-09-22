\name{stackRasterList}
\alias{stackRasterList}
\alias{spdfToBrick}

\title{
Converts a list of rasters, possibly with different projections and resolutions, to a single raster stack.
}
\description{
This function is intended to be used prior to passing covariates to \link{krigeLgm} in order
for the rasters for all covariates to have the same projection and same resolution.  

}
\usage{
stackRasterList(x, template = x[[1]], method = "near", mc.cores=NULL)
spdfToBrick(x, 
    template,
    logSumExpected=FALSE,
    pattern = '^expected_[[:digit:]]+$'
)
}
\arguments{
  \item{x}{
A list of \code{SpatRaster} or \code{SpatVector}s for 
\code{stackRasterList} and \code{spdfToBrick} respectively
}
\item{template}{
A raster whose projection and resolution all other rasters will be aligned with.
}
\item{method}{
The method to use, either "near", or "bilinear".  Can be a vector of the same length as x to specify different methods for each raster.  If \code{method} has names 
which correspond to the names of \code{x}, the names will be used instead of the order to assign methods to rasters.
}
\item{mc.cores}{If non-null, \code{\link[parallel:parallel-package]{mclapply}} is used with this argument specifying the number of cores.}
\item{logSumExpected}{return the log of the sum of offsets}
\item{pattern}{expression to identify layers to rasterize in \code{x}}
}

\value{
A raster brick, with one layer for each variable.
}


\examples{

myCrs = crs("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
x = list(a=rast(matrix(1:9, 3, 3), extent=ext(0,1,0,1), 
       crs=myCrs),
	b=rast(matrix(1:25, 5, 5), extent=ext(-1, 2, -1, 2),
	   crs=myCrs)
	)
mystack = stackRasterList(x)
mystack



mylist = list(
	a=rast(matrix(1:36, 6, 6,byrow=TRUE), extent=ext(0,1000,0,1000), 
       crs=myCrs),
	b=rast(matrix(1:144, 12, 12), extent=ext(-200, 200, -200, 200),
	   crs=myCrs),
	c=rast(matrix(1:100, 10, 10), extent=ext(-5000,5000,-5000,5000), 
       crs=myCrs)
	)
	
	mystack = stackRasterList(mylist, mc.cores=1)
	mystack

plot(mystack[["b"]], main="stack b")
plot(mystack[['a']],add=TRUE,col=grey(seq(0,1,len=12)),alpha=0.8,legend=FALSE)

}

