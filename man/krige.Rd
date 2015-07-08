
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
A SpatialPointsDataFrame containing  the data to be interpolated
}
\item{grid}{Either a \code{\link[raster]{raster}}, or a single integer giving the
number of cells in the X direction which predictions will be made on.  If the later
the predictions will be a raster of square cells covering the bounding box of \code{data}.}
 \item{covariates}{
The spatial covariates used in prediction, either a \code{\link[raster]{raster}} stack or list of rasters.
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
\item{mc.cores}{passed to \code{\link[parallel]{mcmapply}} if greater than 1.}
}
\details{
Given the model parameters and observed data, conditional means and variances of the spatial random field
are computed.
}
\value{
A raster stack is 
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
swissRain$lograin = log(swissRain$rain)
swissRain[[names(swissAltitude)]] = extract(swissAltitude, swissRain)

if(interactive()  | Sys.info()['user'] =='patrick') {
swissFit = likfitLgm(data=swissRain, 
			formula=lograin~ CHE_alt,
			param=c(range=46500, nugget=0.05,shape=1,  
					anisoAngleDegrees=35, anisoRatio=12),
			paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio")
)
myTrend = swissFit$model$formula
myParams = swissFit$param
dput(myParams)
# will give the following
} else {
myParams=structure(c(0.0951677553339405, 0.773131120360615, 49376.366554348, 
1, 11.6731542330799, 0.649923845971017, 1, 2.26103865436198, 
0.000146945159199578, 37.2378933790499), .Names = c("nugget", 
"variance", "range", "shape", "anisoRatio", "anisoAngleRadians", 
"boxcox", "(Intercept)", "CHE_alt", "anisoAngleDegrees"))
myTrend =lograin~ CHE_alt
}

# make sure krige can cope with missing values!
swissAltitude[1:50,1:50] = NA
swissKrige = krigeLgm(data=swissRain, 
	formula = myTrend,
	covariates = swissAltitude,  
	param=myParams,
	grid = 40, expPred=TRUE)



plot(swissKrige[["predict"]], main="predicted rain")
plot(swissBorder, add=TRUE)

# now with box cox and provide a raster for prediction, no covariates

if(interactive()  | Sys.info()['user'] =='patrick') {
swissFit2 = likfitLgm(data=swissRain, 
			formula=rain~1,
			param=c(range=52000, nugget=0.1,shape=1, boxcox=0.5,
					anisoAngleDegrees=35, anisoRatio=8),
			paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"),
			parscale = c(range=5000,nugget=0.01, 
				anisoRatio=1,anisoAngleDegrees=5)
)
myTrend2 = swissFit2$trend
myParams2 = swissFit2$param
dput(myParams2)
} else {
myParams2=structure(c(0.865530531647866, 8.76993204385615, 54143.5826959284, 
1, 7.36559089705556, 0.647158492167979, 0.5, 5.16254700135706, 
37.0794502772753), .Names = c("nugget", "variance", "range", 
"shape", "anisoRatio", "anisoAngleRadians", "boxcox", "(Intercept)", 
"anisoAngleDegrees"))
myTrend2=rain~1
}

swissRaster = raster(extent(swissBorder), nrows=25, ncols=40)


swissKrige2 = krigeLgm(data=swissRain, formula = myTrend2,
	  param=myParams2,
	grid = swissRaster)





plot(swissKrige2[["predict"]], main="predicted rain with box-cox")
plot(swissBorder, add=TRUE)




}
