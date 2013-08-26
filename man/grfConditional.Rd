\name{grfConditional}
\alias{grfConditional}

\title{
Conditional simulation of a Gaussian Random Field
}
\description{
Simulates a Gaussian random field conditional on observed values.
}
\usage{
grfConditional(data,  trend, param, locations, Nsim, 
	covariates,   fun, nugget.in.prediction=TRUE) 
}
 
\arguments{
\item{data}{ A SpatialPointsDataFrame object containing the data to be conditioned upon.
}
\item{trend}{Either a model formula, or a data frame of linear covariates.}

\item{param}{A vector of named model parameters, as produced by \code{\link{likfitLgm}}
}

\item{locations}{Either a \code{\link[raster]{raster}}, or a single integer giving the
number of cells in the X direction which the field  will simulated on.  If the later
the predictions will be a raster of square cells covering the bounding box of \code{locations}.}
\item{Nsim}{Number of samples to simulate.}
  \item{covariates}{
The spatial covariates used in prediction, either a \code{\link[raster]{raster}} stack or list of rasters.
}
 \item{fun}{A function applied to each realised surface and only the result saved.}
\item{nugget.in.prediction}{If TRUE, predict new observations by adding the nugget effect.  Otherwise
predict fitted values.  Only relevant for Box-Cox or log transformed data.}
}

\value{
\item{random}{the predictions and conditional variance of the random effects, on the same
raster as \code{locations}}
}


\seealso{
\code{\link{krige}}
}
\examples{
data(swissRain)
swissRain$elevation = extract(swissAltitude, swissRain)
swissRain$sqrtrain = sqrt(swissRain$rain)

swissFit =  likfitLgm(swissRain, trend=sqrtrain ~ elevation, 
		param=c(range=51000, nugget=0.1,rough=1,  
					aniso.angle.degrees=35, aniso.ratio=7),
		paramToEstimate = c("range","nugget", 
				"aniso.angle.degrees", "aniso.ratio"))


swissRaster = raster(extent(swissBorder), ncols=20, nrows=20, 
	crs=swissRain@proj4string)	

swissSim = grfConditional(data=swissRain, param=swissFit$param,
trend=swissFit$trend, 
locations=swissRaster, Nsim=2, 
covariates=swissAltitude, nugget.in.prediction=TRUE) 

plot(mask(swissSim[[2]], swissBorder)^2) 
plot(swissBorder, add=TRUE)

# now don't include nugget effect
swissSimNoNugget = grfConditional(data=swissRain, param=swissFit$param,
trend=swissFit$trend, 
locations=swissRaster, Nsim=2, 
covariates=swissAltitude, nugget.in.prediction=FALSE) 

plot(mask(swissSimNoNugget[[2]], swissBorder)^2) 
plot(swissBorder, add=TRUE)

# proportion the realisation above 200
swissSim2 = grfConditional(
data=swissRain, param=swissFit$param,
trend=swissFit$trend, 
locations=swissRaster, Nsim=5, 
covariates=swissAltitude, nugget.in.prediction=FALSE,
fun = function(x) mean(x^2>200)  )


unlist(swissSim2) 


}
