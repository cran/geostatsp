\name{lgm-methods}
\docType{methods}
\alias{lgm}
\alias{lgm-methods}
\alias{lgm,ANY,Raster-method}
\alias{lgm,ANY,Spatial-method}

\title{
 Linear Geostatistical Models
}
\description{
Calculate MLE's of model parameters and perform spatial prediction. }
 
 

\usage{
\S4method{lgm}{ANY,Spatial}(formula, data,covariates=NULL, 
		shape=1, boxcox=1, nugget = 0, 
		newdata, 
		expPred=FALSE, nuggetInPrediction=TRUE,
		reml=TRUE,mc.cores=1,
		aniso=FALSE,
		fixShape=TRUE,
		fixBoxcox=TRUE,
		fixNugget = FALSE,
		...) 
\S4method{lgm}{ANY,Raster}(formula, data, 	
covariates=NULL,
                  shape=1,boxcox=1,nugget=0,
				  newdata=data,
				  expPred=FALSE, nuggetInPrediction=TRUE,
				  reml = TRUE,
				  mc.cores=1,		
				  oneminusar=seq(0.01, 0.4,len=4), 
                  range=NULL, 
                  ...) 

}
 
\arguments{
\item{formula}{A model formula for the fixed effects, or a character string specifying the response variable.}
 \item{data}{
A \code{SpatialPointsDataFrame} or \code{Raster} layer, brick or stack containing the locations and observations, and possibly covariates.
}
\item{covariates}{
The spatial covariates used in prediction, either a \code{\link[raster]{raster}} stack or list of rasters.
Covariates in \code{formula} but not in \code{data} will be extracted from \code{covariates}.
}
\item{shape}{Order of the Matern correlation}
\item{boxcox}{Box-Cox transformation parameter (or vector of parameters), set to 1 for no transformation.}
\item{nugget}{Value for the nugget effect (observation error) variance, or vector of such values.}
\item{newdata}{Either a \code{\link[raster]{raster}}, or a single integer giving the
number of cells in the X direction which predictions will be made on.  If the later
the predictions will be a raster of square cells covering the bounding box of \code{data}.}
 \item{expPred}{
Should the predictions be exponentiated, defaults to \code{FALSE}.
}
\item{nuggetInPrediction}{If \code{TRUE}, predict new observations by adding the 
nugget effect.  The prediction variances will be adjusted accordingly, and the predictions 
on the natural scale for logged or Box Cox transformed data will be affected. 
Otherwise predict fitted values. 
}
\item{reml}{If \code{TRUE} (the default), use restricted maximum likelihood.}
\item{mc.cores}{If \code{mc.cores>1}, this argument is passed to  \code{\link[parallel]{mcmapply}} and computations are 
done in parallel where possible. }
\item{aniso}{Set to \code{TRUE} to use geometric anisotropy.}
\item{fixShape}{Set to \code{FALSE} to estimate the Matern order}
\item{fixBoxcox}{Set to \code{FALSE} to estimate the Box-Cox parameter.}
\item{fixNugget}{Set to \code{FALSE} to estimate the nugget effect parameter.}
\item{oneminusar}{Vector of GMRF autoregressive parameters over which the likelihood is computed.}
\item{range}{Vector of range parameters over which the likelihood is computed, ignored if \code{oneminusar} is non-NULL}
\item{...}{Additional arguments passed to \code{\link{likfitLgm}}.  Starting values can be
specified with a vector \code{param} of named elements}
}
\details{
When \code{data} is a \code{SpatialPointsDataFrame}, parameters are estimated using \code{\link[stats]{optim}} to maximize
 the
log-likelihood function computed by  \code{\link{likfitLgm}} and spatial prediction accomplished with \code{\link{krigeLgm}}.

With \code{data} being a \code{Raster} object, a Markov Random Field approximation to the Matern is used (experimental).  Parameters to 
be estimated should be provided as vectors of possible values, with optimization only considering the parameter values supplied.
}
\value{
A list is returned which includes a \code{RasterStack} named \code{predict} having layers:
\item{fixed}{Estimated means from the fixed effects portion of the model}
\item{random}{Predicted random effect}
\item{krigeSd}{Conditional standard deviation of predicted random effect (on the transformed scale if applicable)}
\item{predict}{Prediction of the response, sum of predicted fixed and random effects.  
For Box-Cox or log-transformed data on the natural (untransformed) scale.  
}
\item{predict.log}{If \code{exp.pred=TRUE}, the prediction of the logged process.}
\item{predict.boxcox}{If a box cox transformation was used, the prediction of the process on the transformed scale.}

In addition, the element \code{summery} contains a table of parameter estimates and confidence intervals.  \code{optim} contains the 
output from the call to the  \code{\link[stats]{optim}} function.
}

\seealso{
\code{\link{likfitLgm}}, \code{\link{krigeLgm}}
}
\examples{

data("swissRain")
 
swissRes =  lgm( formula="rain", 
	data=swissRain, newdata=20,
	covariates=swissAltitude, boxcox=0.5, fixBoxcox=TRUE, 
	shape=1, fixShape=TRUE,
	aniso=FALSE, nugget=0, fixNugget=TRUE,
	nuggetInPrediction=FALSE
)

swissRes$summary

plot(swissRes$predict[["predict"]], main="predicted rain") 
plot(swissBorder, add=TRUE)

data("swissRainR")
swissRainR2 = brick(swissRainR[['alt']], 
	 sqrt(swissRainR[['prec1']]))

if(interactive() | Sys.info()['user'] =='patrick') {	 
swissResR =  lgm( formula=layer ~ alt, 
	data=swissRainR2, shape=1,
	oneminusar=seq(0.06, 0.09, len=2),
	nugget =  seq(0,0.0001,len=2)
	)

swissResR$summary[c('oneminusar','range','propNugget','(Intercept).betaHat', 'x.betaHat'),]
}


\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))

library('RColorBrewer')
par(mar=c(0,0,0,3))
plot(gadm) 
plot(mask(projectRaster(
	swissRes$predict[["predict"]],  crs=gadm@proj4string),gadm), 
add=T,alpha=0.6, col=brewer.pal(9, "Blues"))
plot(gadm, add=TRUE)

}




}

 