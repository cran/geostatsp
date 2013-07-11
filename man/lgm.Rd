\name{lgm}
\alias{lgm}
 

\title{
 Linear Geostatistical Models
}
\description{
Fits a   linear geostatistical model using Maximum Likelihood Estimation }
\usage{
lgm(data,  cells, covariates=NULL, formula=NULL,
maternRoughness=1,fixMaternRoughness=TRUE, 
aniso=FALSE, boxcox=1, fixBoxcox=TRUE,
 nugget = 0, fixNugget = FALSE,...)
}
 
\arguments{
\item{data}{
An object of class \code{ SpatialPointsDataFrame} containing the data.
}
\item{cells}{Either an integer giving the number of cells in the x direction, or a raster object which 
will be used for the spatial random effect.  If the cells in the raster are not square, the resolution in the y direction 
will be adjusted to make it so.}
\item{covariates}{Either a single raster, a list of rasters or a raster stack containing covariate values used when 
making spatial predictions.  Names of the raster layers or list elements correspond to names in the formula.  If
a covariate is missing from the data object it will be extracted from the rasters.  Defaults to \code{NULL} for an 
intercept-only model.
}
\item{formula}{Fixed effects formula, 
defaults to a linear combination of each of the layers in the \code{covariates} object.   
 The response variable defaults to the first variable in 
 the \code{data} object, and
 \code{formula} can be an integer or character string specifying the response variable. 
}
\item{maternRoughness}{Order of the Matern correlation}
\item{fixMaternRoughness}{Set to \code{FALSE} to estimate the Matern order}
\item{aniso}{Set to \code{TRUE} to use geometric anisotropy.}
\item{boxcox}{Box-Cox transformation parameter, set to 1 for no transformation.}
\item{fixBoxcox}{Set to \code{FALSE} to estimate the Box-Cox parameter.}
\item{nugget}{Value for the nugget effect (observation error}
\item{fixNugget}{Set to \code{FALSE} to estimate the nugget effect parameter.}
\item{...}{Additional options passed to  \code{\link{likfit}}}
}
\details{
Calls  \code{\link{likfit}} and  \code{\link{krige}}
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
raster as \code{cells}}
Additionally, returned are
\item{parameters}{A data frame with parameter estimates, and standard errors of the fixed
effects parameters}
\item{likfit}{results from \code{\link{likfit}}}
}

\seealso{
\code{\link{likfit}}, \code{\link{krige}}
}
\examples{

data(swissRain)
 
swissRes =  lgm(swissRain, cells=30, formula="rain",
	covariates=swissAltitude, boxcox=0.5, fixBoxcox=TRUE, 
	aniso=TRUE, 
	# starting values, arguments as in geoR::likfit
	ini.cov.pars=c(70, 20000), psiA=0.65, psiR=8
	)

swissRes$parameters

plot(swissRes$predict[["predict"]]) 
\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))
plot(spTransform(gadm, swissRain@proj4string), add=TRUE)
}

}


 	