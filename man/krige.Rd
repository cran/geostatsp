\name{krige}
\alias{krige}

\title{
Spatial prediction, or Kriging
}
\description{
Perform spatial prediction, producing a raster of predictions and conditional standard deviations.
}
\usage{
krige(data, trend, coordinates = data, param, locations, 
    covariates = list(), expPred = FALSE, nugget.in.prediction = TRUE,
    conditionalVarianceMatrix=FALSE)
}
 
\arguments{
  \item{data}{
A SpatialPointsDataFrame containing  the data to be interpolated, or a vector of observed data
}
\item{trend}{Either a model formula, or a data frame of linear covariates.}

\item{coordinates}{A SpatialPoints object giving the locations of the observed data.  
If \code{data} is a SpatialPointsDataFrame coordinates is not necessary.}
\item{covariates}{
The spatial covariates used in prediction, either a \code{\link[raster]{raster}} stack or list of rasters.
}
\item{param}{A vector of named model parameters, as produced by \code{\link{likfitLgm}}
}
\item{locations}{Either a \code{\link[raster]{raster}}, or a single integer giving the
number of cells in the X direction which predictions will be made on.  If the later
the predictions will be a raster of square cells covering the bounding box of \code{locations}.}
  \item{expPred}{
Should the predictions be exponentiated, defaults to \code{FALSE}.
}
\item{nugget.in.prediction}{If \code{TRUE}, predict new observations by adding the 
nugget effect.  The prediction variances will be adjusted accordingly, and the predictions 
on the natural scale for logged or Box Cox transformed data will be affected. 
Otherwise predict fitted values. 
}
\item{conditionalVarianceMatrix}{Return the full conditional variance matrix.}
}
\details{
Given the model parameters and observed data, conditional means and variances of the spatial random field
are computed.
}
\value{
If the prediction locations are identical for the fixed effects and random effects, a single raster stack is 
returned with the following layers:
\item{fixed}{Estimated means from the fixed effects portion of the model}
\item{random}{Predicted random effect}
\item{krige.var}{Conditional variance of predicted random effect (on the transformed scale if applicable)}
\item{predict}{Prediction of the response, sum of fixed and random effects.  If exp.pred is TRUE, gives predictions
on the exponentated scale, and half of krige.var is added prior to exponentiating}
\item{predict.log}{If exp.pred=TRUE, the prediction of the logged process.}
If the prediction locations are different for fixed and random effects (typically coarser for the random 
effects), a list with two raster stacks is returned.
\item{predict.boxcox}{If a box cox transformation was used, the prediction of the process on the transformed scale.}
If the prediction locations are different for fixed and random effects (typically coarser for the random 
effects), a list with two raster stacks is returned.
\item{prediction}{A raster stack as above, though the random effect prediction is resampled to the 
same locations as the fixed effects.}
\item{random}{the predictions and conditional variance of the random effects, on the same
raster as \code{locations}}


}


\seealso{
\code{\link{lgm}}
}
\examples{
data(swissRain)
swissRain$lograin = log(swissRain$rain)
swissRain$elevation = extract(swissAltitude, swissRain)

swissFit = likfitLgm(data=swissRain, 
			trend=lograin~ elevation,
			param=c(range=46500, nugget=0.05,rough=1,  
					aniso.angle.degrees=35, aniso.ratio=12),
			paramToEstimate = c("range","nugget", 
				"aniso.angle.degrees", "aniso.ratio"),
			parscale = c(range=5000,nugget=0.01, 
				aniso.ratio=1,aniso.angle.degrees=5)
)

# make sure krige can cope with missing values!
swissAltitude[1:50,1:50] = NA
swissKrige = krige(data=swissRain, trend = swissFit$trend,
	covariates = swissAltitude,  
	param=swissFit$param,
	locations = 40)

plot(swissKrige[["predict"]])
plot(swissBorder, add=TRUE)

# now with box cox and provide a raster for prediction

swissFit2 = likfitLgm(data=swissRain, 
			trend=rain~1,
			param=c(range=52000, nugget=0.1,rough=1, boxcox=0.5,
					aniso.angle.degrees=35, aniso.ratio=8),
			paramToEstimate = c("range","nugget", 
				"aniso.angle.degrees", "aniso.ratio"),
			parscale = c(range=5000,nugget=0.01, 
				aniso.ratio=1,aniso.angle.degrees=5)
)
swissRaster = raster(extent(swissBorder), nrows=25, ncols=40)


swissKrige2 = krige(data=swissRain, trend = swissFit2$trend,
	  param=swissFit2$param,
	locations = swissRaster, expPred=TRUE)

plot(swissKrige2[["predict"]])
plot(swissBorder, add=TRUE)




}
