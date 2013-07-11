\name{simLgcp}
\alias{simLgcp}
\alias{simPoissonPP}


\title{Simulate a log-Gaussian Cox process}
\description{
Give covariates and model parameters, simulates a log-Gaussian Cox process
}
\usage{
	simLgcp(param, covariates=NULL, betas=NULL, 
	  rasterTemplate=covariates[[1]], model="whittle", ...)
	simPoissonPP(intensity)
}

\arguments{
\item{param}{Model parameters for the latent Gaussian surface passed to \code{\link[RandomFields]{GaussRF}}.  For a Matern model, \code{params} should be a vector with elements \code{mean},
\code{variance}, \code{scale}, \code{alpha}}.
\item{covariates}{Either a raster stack or list of rasters and SpatialPolygonsDataFrames (with the latter having only a single data column).}
\item{betas}{Coefficients for the covariates}
\item{rasterTemplate}{Raster on which the latent surface is simulated,
defaults to the first covariate.}
\item{model}{model passed to \code{\link[RandomFields]{GaussRF}}}, defaults to \code{whittle} for a Matern correlation.
\item{...}{ additional arguments, see \code{\link[RandomFields]{GaussRF}} in the RandomFields package.}
\item{intensity}{Raster of the intensity of a Poisson point process.}
}

\value{
A list with a \code{events} element containing the event locations and a \code{raster} element 
containing a raster stack of the covariates, spatial random effect, and intensity.
}

\examples{
mymodel = c(mean=-0.5, variance=1, 
				scale=2, alpha=2)

myraster = raster(nrows=20,ncols=20,xmn=0,xmx=10,ymn=0,ymx=10)

myCovariate = list(a=raster(matrix(1:10,10,10), 0,10,0,10),
b=raster(matrix(1:10,10,10,byrow=TRUE), 0,10,0,10)) 

myLgcp=simLgcp(mymodel, myCovariate, betas=c(a=-0.1, b=0.25), 
	rasterTemplate=myraster)

plot(myLgcp$raster[["intensity"]])
points(myLgcp$events)

myIntensity = exp(-1+0.2*myCovariate[["a"]])
myPoissonPP = simPoissonPP(myIntensity)
plot(myIntensity)
points(myPoissonPP)


}

