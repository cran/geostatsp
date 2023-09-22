
\name{krigeLgm}
\alias{krigeLgm}

\title{
Spatial prediction, or Kriging
}
\description{
Perform spatial prediction, producing a raster of predictions and conditional standard deviations.
}
\usage{
krigeLgm(formula, data, grid,  covariates = NULL,
	param,  
    expPred = FALSE, nuggetInPrediction = TRUE,
    mc.cores=getOption("mc.cores", 1L))
}
 
\arguments{
\item{formula}{Either a model formula, or a data frame of linear covariates.}
  \item{data}{
A \code{SpatVector}  containing  the data to be interpolated
}
\item{grid}{Either a \code{\link[terra]{SpatRaster}} , or a single integer giving the
number of cells in the X direction which predictions will be made on.  If the later
the predictions will be a raster of square cells covering the bounding box of \code{data}.}
 \item{covariates}{
The spatial covariates used in prediction, either a \code{\link[terra]{SpatRaster}} stack or list of rasters.
}
\item{param}{A vector of named model parameters, as produced by \code{\link{likfitLgm}}
}
 \item{expPred}{
Should the predictions be exponentiated, defaults to \code{FALSE}.
}
\item{nuggetInPrediction}{If \code{TRUE}, predict new observations by adding the 
nugget effect.  The prediction variances will be adjusted accordingly, and the predictions 
on the natural scale for logged or Box Cox transformed data will be affected. 
Otherwise predict fitted values. 
}
\item{mc.cores}{passed to \code{\link[parallel:parallel-package]{mclapply}} if greater than 1.}
}
\details{
Given the model parameters and observed data, conditional means and variances of the spatial random field
are computed.
}
\value{
A raster  is 
returned with the following layers:
\item{fixed}{Estimated means from the fixed effects portion of the model}
\item{random}{Predicted random effect}
\item{krige.var}{Conditional variance of predicted random effect (on the transformed scale if applicable)}
\item{predict}{Prediction of the response, sum of fixed and random effects.  If exp.pred is TRUE, gives predictions
on the exponentiated scale, and half of krige.var is added prior to exponentiating}
\item{predict.log}{If exp.pred=TRUE, the prediction of the logged process.}
\item{predict.boxcox}{If a box cox transformation was used, the prediction of the process on the transformed scale.}
If the prediction locations are different for fixed and random effects (typically coarser for the random 
effects), a list with two raster stacks is returned.
\item{prediction}{A raster stack as above, though the random effect prediction is resampled to the 
same locations as the fixed effects.}
\item{random}{the predictions and conditional variance of the random effects, on the same
raster as \code{newdata}}

}


\seealso{
\code{\link{lgm}}
}
\examples{
data('swissRain')
swissAltitude = unwrap(swissAltitude)
swissRain = unwrap(swissRain)
swissRain$lograin = log(swissRain$rain)
swissRain[[names(swissAltitude)]] = extract(swissAltitude, swissRain, ID=FALSE)

swissFit = likfitLgm(data=swissRain, 
			formula=lograin~ CHE_alt,
			param=c(range=46500, nugget=0.05,shape=1,  
					anisoAngleDegrees=35, anisoRatio=12),
			paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio")
)
myTrend = swissFit$model$formula
myParams = swissFit$param


swissBorder = unwrap(swissBorder)

swissKrige = krigeLgm(
	data=swissRain, 
	formula = myTrend,
	covariates = swissAltitude,  
	param=myParams,
	grid = squareRaster(swissBorder, 40), expPred=TRUE)

plot(swissKrige[["predict"]], main="predicted rain")
plot(swissBorder, add=TRUE)




}
