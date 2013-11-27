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
	GaussRF(x, param=c(variance=1, range=1, shape=1), ...)
	\method{GaussRF}{Raster}(x,param=c(variance=1, range=1, shape=1),...)
	\method{GaussRF}{SpatialPoints}(x,param=c(variance=1, range=1, shape=1),...)
	\method{GaussRF}{SpatialPointsDataFrame}(x,param=c(variance=1, range=1, shape=1),...)
}

\arguments{
	\item{x}{An object of class \code{\link[raster]{raster}}, SpatialPoints, or SpatialPointsDataFrame,
	or an object of a class acceptable to the  \code{\link[RandomFields]{GaussRF}} .}
	\item{param}{A vector of named model parameters with, at a minimum names
	\code{variance},  \code{range} and \code{shape}.  For Geometric Anisotropy add 
	\code{anisoRatio} and either \code{anisoAngleDegrees} or  
		\code{anisoAngleRadians}
	}
	\item{...}{ additional arguments, see \code{\link[RandomFields]{GaussRF}} in the RandomFields package.}
}

\value{
	If \code{x} is of class  \code{\link[raster]{raster}}, a random field is simulated on a grid and returned as a raster. 
	 If \code{x} is of class Spatial
	Points or SpatialPointsDataFrame, a random field is simulated at the points specified and a vector is returned.

}
\details{
The formula for the Matern correlation function is
\deqn{M(x) =  \frac{variance}{\Gamma(shape)} 2^{shape-1} 
x (\sqrt{8 shape}/ range)^{shape}  
besselK(x \sqrt{8 shape}/ range, shape) }{%
             M(x) =   (variance / Gamma(shape)) 2^(shape-1) 
x [sqrt(8 shape)/ range]^shape  
besselK[x sqrt(8 shape)/ range, shape] 
}

The \code{range} parameter is sqrt(8*shape)*phi (for phi in the geoR package), 
sqrt(8*shape)*scale (for 'scale' as in the whittle model in RandomFields), and
2*scale (for 'scale' in the matern model in RandomFields).

If an argument \code{model} is specified, parameters are passed directly to the RandomFields package and 
	the parametrizations from that package are used.
}

\examples{
mymodel = c(variance=1, range=2, shape=2, anisoRatio=0.6, anisoAngleDegrees=30)

myraster = raster(nrows=40,ncols=40,xmn=-20,xmx=20,ymn=-15,ymx=15)

myGRF=GaussRF(myraster, param=mymodel)

plot(myGRF, main="a GRF")

mypoints = SpatialPointsDataFrame(cbind(runif(10), runif(10)),data=data.frame(id=1:10))

mypoints$rf= GaussRF(mypoints, param=mymodel) 

spplot(mypoints, "rf", col.regions=terrain.colors(5), main="GRF at points")


x = cbind(sample(1:50), sample(1:50))
temp=GaussRF(x,  
		param=c(variance=1,	range=6, shape=2))
quantile(temp)

# passing parameters straight to RandomFields
fromRF = GaussRF(myraster, param=c(mean=0.5, variance=1, nugget=0, 
				scale=4, alpha=2), model="matern")
plot(fromRF, main="from RandomFields")

}
