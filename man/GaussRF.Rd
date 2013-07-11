\name{GaussRF}
\alias{GaussRF}
\alias{GaussRF.default}
\alias{GaussRF.Raster}
\alias{GaussRF.SpatialPoints}
\alias{GaussRF.SpatialPointsDataFrame}

\title{Simulate Gaussian Random Fields}
\description{
	This function is a wrapper for the \code{\link[RandomFields]{GaussRF}} 
function in the \link[RandomFields]{RandomFields} package, where the arguments and output
are Raster or SpatialPoints objects.
}
\usage{
	GaussRF(x,  ...)
	\method{GaussRF}{Raster}(x,...)
	\method{GaussRF}{SpatialPoints}(x,...)
	\method{GaussRF}{SpatialPointsDataFrame}(x,...)
}

\arguments{
	\item{x}{An object of class \code{\link[raster]{raster}}, SpatialPoints, or SpatialPointsDataFrame,
	or an object of a class acceptable to the  \code{\link[RandomFields]{GaussRF}} .}
\item{...}{ additional arguments, see \code{\link[RandomFields]{GaussRF}} in the RandomFields package.}
}

\value{
	If \code{x} is of class  \code{\link[raster]{raster}}, a random field is simulated on a grid and returned as a raster. 
	 If \code{x} is of class Spatial
	Points or SpatialPointsDataFrame, a random field is simulated at the points specified and a vector is returned.
}

\examples{
mymodel = c(mean=0, variance=1, nugget=0, 
				scale=2, alpha=2)

myraster = raster(nrows=11,ncols=11,xmn=0,xmx=10,ymn=0,ymx=10)

myGRF=GaussRF(myraster, model="whittle", 
		param=mymodel)

plot(myGRF)

mypoints = SpatialPointsDataFrame(cbind(runif(10), runif(10)),data=data.frame(id=1:10))

mypoints$rf= GaussRF(mypoints, model="whittle", 
		param=mymodel) 

spplot(mypoints, "rf")


x = cbind(1:4, 1:4)
temp=GaussRF(x, model="whittle", 
		param=c(mean=0, variance=1, nugget=0, 
				scale=2, alpha=2))


}
