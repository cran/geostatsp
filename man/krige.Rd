\name{krige}
\alias{krige}

\title{
Spatial prediction, or Kriging
}
\description{
This function uses the  \code{\link[geoR]{krige.control}} and  \code{\link[geoR]{krige.conv}} functions
in geoR to perform spatial prediction using a fitted model from \code{\link{likfit}}.
}
\usage{
krige(obj.model, geodata, locations, covariates, locations.mean = locations,
   factor.info = NULL, 
	exp.pred = FALSE, rasterMethod = c("ngb", "bilinear"),
	nugget.in.prediction=TRUE, ...)
}
 
\arguments{
  \item{obj.model}{
A fitted model produced by \code{\link{likfit}}
}
  \item{geodata}{
A SpatialPointsDataFrame containing  the data used with \code{\link{likfit}}.
}
  \item{locations}{
A \code{\link[raster]{raster}} on which the spatial random effects will be predicted, or an integer specifying
the number of cells in the x dimension for a raster with square cells covering \code{geodata}.
}
  \item{covariates}{
The spatial covariates used in prediction, either a \code{\link[raster]{raster}} stack or list of rasters.
}
  \item{locations.mean}{
A \code{\link[raster]{raster}} on which the fixed effects component of the model is predicted.  Defaults
to the same raster as used for the random effects.  Ignored if \code{covariates} is a raster stack.
}
  \item{factor.info}{
A list of data frames, with first column being levels and second column being labels for covariates which are 
factors.  Useful if the rasters containing categorical covariates are integers corresponding to different 
categories from those used in \code{geodata}.
}
  \item{exp.pred}{
Should the predictions be exponentiated, defaults to false.
}
  \item{rasterMethod}{
Method for reprojecting or resampling rasters, either "bilinear" or "ngb". Any covariates which are factors
will be reprojected with ngb even if bilinear is specified.
}
\item{nugget.in.prediction}{If TRUE, predict new observations by adding the nugget effect.  Otherwise
predict fitted values.  Only relevant for Box-Cox or log transformed data.
}
\item{...}{Additional arguments, currently none
}
}
\details{
This function is essentially a wrapper for \code{\link[geoR]{krige.conv}}.
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
\code{\link[geoR]{krige.conv}}
}
\examples{
data(swissRain)
swissRain$lograin = log(swissRain$rain)
swissRain$elevation = extract(swissAltitude, swissRain)
swissFit =  likfit(swissRain, lograin ~ elevation, 
	cov.model="matern", kappa = 2, fix.kappa=TRUE)
swissRaster = raster(extent(swissRain@bbox), ncols=61, nrows=91, 
  crs=swissRain@proj4string)	
swissKrige = krige(swissFit, swissRain, swissRaster, swissAltitude)


swissFit =  likfit(swissRain, lograin ~ elevation, 
	cov.model="matern", kappa = 2, fix.kappa=TRUE,
	lambda=0.5, fix.lambda=TRUE,
	ini.cov.pars= c(0.6, 7000))
swissKrige = krige(obj.model=swissFit, geodata=swissRain, 
	locations=swissRaster, covariates=swissAltitude)


plot(swissKrige[["predict"]])
plot(swissBorder, add=TRUE)

 

}

 

