
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

\dontrun{
swissFit = likfitLgm(data=swissRain, 
			formula=lograin~ CHE_alt,
			param=c(range=46500, nugget=0.05,shape=1,  
					anisoAngleDegrees=35, anisoRatio=12),
			paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"),
			parscale = c(range=5000,nugget=0.01, 
				anisoRatio=1,anisoAngleDegrees=5)
)
myTrend = swissFit$formula
myParams = swissFit$param
dput(myParams)
# will give the following
}
myParams=structure(c(2.26621791403682, 0.00014474608227038, 48928.027588747, 
0.0973034174963395, 37.2883307364709, 11.5429109094065, 1, 0.736723861095649
), .Names = c("(Intercept)", "CHE_alt", "range", "nugget", "anisoAngleDegrees", 
"anisoRatio", "shape", "variance"))
myTrend =lograin~ CHE_alt


# make sure krige can cope with missing values!
swissAltitude[1:50,1:50] = NA
swissKrige = krigeLgm(data=swissRain, formula = myTrend,
	covariates = swissAltitude,  
	param=myParams,
	grid = 40, expPred=TRUE)



plot(swissKrige[["predict"]], main="predicted rain")
plot(swissBorder, add=TRUE)

# now with box cox and provide a raster for prediction, no covariates
\dontrun{
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
}
myParams2=structure(c(5.1259476451533, 57230.2319794406, 0.878030061109591, 
37.1176523248768, 7.37809931787162, 1, 0.5, 9.3502154358886), .Names = c("(Intercept)", 
"range", "nugget", "anisoAngleDegrees", "anisoRatio", "shape", 
"boxcox", "variance"))
myTrend2=rain~1

swissRaster = raster(extent(swissBorder), nrows=25, ncols=40)


swissKrige2 = krigeLgm(data=swissRain, formula = myTrend2,
	  param=myParams2,
	grid = swissRaster)





plot(swissKrige2[["predict"]], main="predicted rain with box-cox")
plot(swissBorder, add=TRUE)




}
