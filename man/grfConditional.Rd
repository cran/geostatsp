\name{grfConditional}
\alias{grfConditional}

\title{
Conditional simulation of a Gaussian Random Field
}
\description{
Simulates a Gaussian random field conditional on observed values.
}
\usage{
grfConditional(data, ycol=1, 
			param, locations, Nsim,
		 	fun=NULL, nuggetInPrediction=TRUE) 
}
\arguments{
\item{data}{ Either a SpatialPoints object of observation locations, or a 
SpatialPointsDataFrame with locations and (residuals from) observed data
 at these locations.
}
\item{ycol}{Either a vector of (residuals from) observed data or name or number of 
the column of \code{data} containing the observed data.}

\item{param}{A vector of named model parameters, as produced by \code{\link{likfitLgm}}
}
\item{locations}{Either a \code{\link[raster]{raster}}, or a single integer giving the
number of cells in the X direction which the field  will simulated on.  If the later
the predictions will be a raster of square cells covering the bounding box of \code{locations}.}
\item{Nsim}{Number of samples to simulate.}
 \item{fun}{A function applied to each realised surface and only the result saved.}
\item{nuggetInPrediction}{If TRUE, predict new observations by adding the nugget effect.  Otherwise
predict fitted values.  Only relevant for Box-Cox or log transformed data.}
}

\value{
A raster of a Gaussian random field simulated conditional on the observed data.  Note the
observed data and  simulated data both have mean zero, so residuals from a fitted model
rather than the data themselves should be supplied.
}


\seealso{
\code{\link{krige}}
}
\examples{
\dontrun{ # it's fairly time consuming
data("swissRain")
swissRain$elevation = extract(swissAltitude, swissRain)
swissRain$sqrtrain = sqrt(swissRain$rain)

# estimate parameters
swissFit =  likfitLgm(swissRain, trend=sqrtrain ~ elevation, 
		param=c(range=51700, nugget=0.11,shape=1,  
				anisoAngleDegrees=37, anisoRatio=7.6),
		paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"))

# simulate from the random effect conditional on
#   the observed data
swissSim = grfConditional(data=swissRain, ycol=swissFit$resid,
		param=swissFit$param, locations=30, 
		Nsim=1)

# plot the simulated random effect
plot(swissSim)
plot(swissBorder, add=TRUE)


# create a small raster of elevation data
swissAltSmall = aggregate(swissAltitude,fact=5)
# calculate the fixed effects portion of the rainfall process
rainMean = swissFit$param["(Intercept)"] +
			swissFit$param["elevation"] * swissAltSmall

# define a function to identify the location of maximum rainfall	
maxRainLocation = function(x) {
	rain =  (rainMean + x)^2
	xyFromCell(rain, which.max(rain))
}

# get a conditional sample of three locations of maximum rainfall
swissRain$resid = swissFit$resid
swissLocation = grfConditional(data=swissRain, 
		ycol="resid",
		param=swissFit$param, locations=swissAltSmall, 
		Nsim=3, fun = maxRainLocation)

# convert result from a list to a matrix
swissLocation = matrix(unlist(swissLocation), ncol=2,byrow=TRUE)
# add the locations to the map
points(swissLocation, pch=1:(dim(swissLocation)[1]),col='red')

}
}
