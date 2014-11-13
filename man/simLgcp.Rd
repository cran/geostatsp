\name{simLgcp}
\alias{simLgcp}
\alias{simPoissonPP}


\title{Simulate a log-Gaussian Cox process}
\description{
Give covariates and model parameters, simulates a log-Gaussian Cox process
}
\usage{
	simLgcp(param, covariates=NULL, betas=NULL, 
	offset=NULL,
	  rasterTemplate=covariates[[1]], n=1, ...)
	simPoissonPP(intensity)
}

\arguments{
\item{param}{A vector of named model parameters with, at a minimum names
 \code{range} and  \code{shape} (see Details), and optionally \code{variance} (defaults to 1).  
	For Geometric Anisotropy add 
	\code{anisoRatio} and either \code{anisoAngleDegrees} or  
		\code{anisoAngleRadians}
	}
	\item{covariates}{Either a raster stack or list of rasters and SpatialPolygonsDataFrames (with the latter having only a single data column).}
\item{betas}{Coefficients for the covariates}
\item{offset}{Vector of character strings corresponding to elements of \code{covariates} which are offsets}
\item{rasterTemplate}{Raster on which the latent surface is simulated,
defaults to the first covariate.}
\item{n}{number of realisations to simulate}
\item{...}{ additional arguments, see \code{\link[RandomFields]{GaussRF}} in the RandomFields package.}
\item{intensity}{Raster of the intensity of a Poisson point process.}
}

\value{
A list with a \code{events} element containing the event locations and a \code{raster} element 
containing a raster stack of the covariates, spatial random effect, and intensity.
}

\examples{
mymodel = c(mean=-0.5, variance=1, 
				range=2, shape=2)

myraster = raster(nrows=15,ncols=20,xmn=0,xmx=10,ymn=0,ymx=7.5)

# some covariates, deliberately with a different resolution than myraster
covA = covB = myoffset = raster(extent(myraster), 10, 10)
values(covA) = as.vector(matrix(1:10, 10, 10))
values(covB) = as.vector(matrix(1:10, 10, 10, byrow=TRUE))
values(myoffset) = round(seq(-1, 1, len=ncell(myoffset)))

myCovariate = list(a=covA, b=covB, offsetFooBar = myoffset)

myLgcp=simLgcp(mymodel, myCovariate, betas=c(a=-0.1, b=0.25), 
	offset='offsetFooBar',
	rasterTemplate=myraster)

plot(myLgcp$raster[["intensity"]], main="lgcp")
points(myLgcp$events)

myIntensity = exp(-1+0.2*myCovariate[["a"]])
myPoissonPP = simPoissonPP(myIntensity)[[1]]
plot(myIntensity, main="Poisson pp")
points(myPoissonPP)


}

