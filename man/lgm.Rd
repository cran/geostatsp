\name{lgm}
\alias{lgm}
 

\title{
 Linear Geostatistical Models
}
\description{
Calculate MLE's of model parameters and perform spatial prediction. }
\usage{
lgm(data, locations, covariates = NULL, formula = NULL, shape = 1, 
    fixShape = TRUE, aniso = FALSE, boxcox = 1, fixBoxcox = TRUE, 
    nugget = 0, fixNugget = FALSE, expPred = FALSE, nuggetInPrediction = TRUE,
    ...) 

}
 
\arguments{
 \item{data}{
A SpatialPointsDataFrame containing  the data to be interpolated, or a vector of observed data
}


\item{covariates}{
The spatial covariates used in prediction, either a \code{\link[raster]{raster}} stack or list of rasters.
Covariates in \code{formula} but not in \code{data} will be extracted from \code{covariates}.
}

\item{locations}{Either a \code{\link[raster]{raster}}, or a single integer giving the
number of cells in the X direction which predictions will be made on.  If the later
the predictions will be a raster of square cells covering the bounding box of \code{locations}.}
\item{formula}{A model formula for the fixed effects.}
\item{shape}{Order of the Matern correlation}
\item{fixShape}{Set to \code{FALSE} to estimate the Matern order}
\item{aniso}{Set to \code{TRUE} to use geometric anisotropy.}
\item{boxcox}{Box-Cox transformation parameter, set to 1 for no transformation.}
\item{fixBoxcox}{Set to \code{FALSE} to estimate the Box-Cox parameter.}
\item{nugget}{Value for the nugget effect (observation error}
\item{fixNugget}{Set to \code{FALSE} to estimate the nugget effect parameter.}

 \item{expPred}{
Should the predictions be exponentiated, defaults to \code{FALSE}.
}
\item{nuggetInPrediction}{If \code{TRUE}, predict new observations by adding the 
nugget effect.  The prediction variances will be adjusted accordingly, and the predictions 
on the natural scale for logged or Box Cox transformed data will be affected. 
Otherwise predict fitted values. 
}
\item{...}{Additional arguments passed to \code{\link{likfitLgm}}.  Starting values can be
specified with a vector \code{param} of named elements}
}
\details{
Calls  \code{\link{likfitLgm}} and  \code{\link{krige}}
}
\value{
A raster stack is 
returned with the following layers:
\item{fixed}{Estimated means from the fixed effects portion of the model}
\item{random}{Predicted random effect}
\item{krigeSd}{Conditional standard deviation of predicted random effect (on the transformed scale if applicable)}
\item{predict}{Prediction of the response, sum of predicted fixed and random effects.  
For Box-Cox or log-transformed data on the natural (untransformed) scale.  
}
\item{predict.log}{If \code{exp.pred=TRUE}, the prediction of the logged process.}
\item{predict.boxcox}{If a box cox transformation was used, the prediction of the process on the transformed scale.}
}

\seealso{
\code{\link{likfitLgm}}, \code{\link{krige}}
}
\examples{

data("swissRain")
 
swissRes =  lgm(swissRain, locations=20, formula="rain",
	covariates=swissAltitude, boxcox=0.5, fixBoxcox=TRUE, 
	shape=1, fixShape=TRUE,
	aniso=FALSE, nugget=0, fixNugget=TRUE,
	nuggetInPrediction=FALSE
	)
	

swissRes$summary

plot(swissRes$predict[["predict"]], main="predicted rain") 
plot(swissBorder, add=TRUE)



\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))

library(RColorBrewer)
par(mar=c(0,0,0,3))
plot(gadm) 
plot(mask(projectRaster(
	swissRes$predict[["predict"]],  crs=gadm@proj4string),gadm), 
add=T,alpha=0.6, col=brewer.pal(9, "Blues"))
plot(gadm, add=TRUE)

}



}

 