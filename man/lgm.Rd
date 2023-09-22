\name{lgm-methods}
\docType{methods}
\alias{lgm}
\alias{lgm-methods}
\alias{lgm,missing,ANY,ANY,ANY-method}
\alias{lgm,numeric,ANY,ANY,ANY-method}
\alias{lgm,character,ANY,ANY,ANY-method}
\alias{lgm,formula,SpatVector,numeric,ANY-method}
\alias{lgm,formula,SpatVector,SpatRaster,missing-method}
\alias{lgm,formula,SpatVector,SpatRaster,list-method}
\alias{lgm,formula,SpatVector,SpatRaster,SpatRaster-method}
\alias{lgm,formula,SpatVector,SpatRaster,data.frame-method}
\alias{lgm,formula,SpatRaster,ANY,ANY-method}
\alias{lgm,formula,data.frame,SpatRaster,data.frame-method}

\title{
 Linear Geostatistical Models
}
\description{
Calculate MLE's of model parameters and perform spatial prediction. }
 
 

\usage{
\S4method{lgm}{missing,ANY,ANY,ANY}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{numeric,ANY,ANY,ANY}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{character,ANY,ANY,ANY}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{formula,SpatVector,numeric,ANY}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{formula,SpatVector,SpatRaster,missing}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{formula,SpatVector,SpatRaster,list}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{formula,SpatVector,SpatRaster,SpatRaster}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{formula,SpatVector,SpatRaster,data.frame}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{formula,SpatRaster,ANY,ANY}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
\S4method{lgm}{formula,data.frame,SpatRaster,data.frame}(
formula, data, grid, covariates, 
buffer=0, shape=1, boxcox=1, nugget = 0, 
expPred=FALSE, nuggetInPrediction=TRUE,
reml=TRUE,mc.cores=1,
aniso=FALSE,
fixShape=TRUE,
fixBoxcox=TRUE,
fixNugget = FALSE,
...) 
}
 
\arguments{
\item{formula}{A model formula for the fixed effects, or a character string specifying the response variable.}
 \item{data}{
A \code{SpatVector} or \code{SpatRaster} layer, brick or stack containing the locations and observations, and possibly covariates.
}
\item{grid}{Either a \code{\link[terra]{SpatRaster}}, or a single integer giving the
number of cells in the X direction which predictions will be made on.  If the later
the predictions will be a raster of square cells covering the bounding box of \code{data}.}
\item{covariates}{
The spatial covariates used in prediction, either a \code{\link[terra]{SpatRaster}} stack or list of rasters.
Covariates in \code{formula} but not in \code{data} will be extracted from \code{covariates}.
}
\item{shape}{Order of the Matern correlation}
\item{boxcox}{Box-Cox transformation parameter (or vector of parameters), set to 1 for no transformation.}
\item{nugget}{Value for the nugget effect (observation error) variance, or vector of such values.}
 \item{expPred}{
Should the predictions be exponentiated, defaults to \code{FALSE}.
}
\item{nuggetInPrediction}{If \code{TRUE}, predict new observations by adding the 
nugget effect.  The prediction variances will be adjusted accordingly, and the predictions 
on the natural scale for logged or Box Cox transformed data will be affected. 
Otherwise predict fitted values. 
}
\item{reml}{If \code{TRUE} (the default), use restricted maximum likelihood.}
\item{mc.cores}{If \code{mc.cores>1}, this argument is passed to  \code{\link[parallel:parallel-package]{mclapply}} and computations are 
done in parallel where possible. }
\item{aniso}{Set to \code{TRUE} to use geometric anisotropy.}
\item{fixShape}{Set to \code{FALSE} to estimate the Matern order}
\item{fixBoxcox}{Set to \code{FALSE} to estimate the Box-Cox parameter.}
\item{fixNugget}{Set to \code{FALSE} to estimate the nugget effect parameter.}
\item{buffer}{Extra distance to add around \code{grid}.}
\item{...}{Additional arguments passed to \code{\link{likfitLgm}}.  Starting values can be
specified with a vector \code{param} of named elements}
}
\details{
When \code{data} is a \code{SpatVector}, parameters are estimated using \code{\link[stats]{optim}} to maximize
 the
log-likelihood function computed by  \code{\link{likfitLgm}} and spatial prediction accomplished with \code{\link{krigeLgm}}.

With \code{data} being a \code{Raster} object, a Markov Random Field approximation to the Matern is used (experimental).  Parameters to 
be estimated should be provided as vectors of possible values, with optimization only considering the parameter values supplied.
}
\value{
A list is returned which includes a \code{SpatRaster} named \code{predict} having layers:
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
swissRain = unwrap(swissRain)
swissAltitude = unwrap(swissAltitude)
swissBorder = unwrap(swissBorder)
 
swissRes =  lgm( formula="rain", 
	data=swissRain[1:60,], grid=20,
	covariates=swissAltitude, boxcox=0.5, fixBoxcox=TRUE, 
	shape=1, fixShape=TRUE,
	aniso=FALSE, nugget=0, fixNugget=FALSE,
	nuggetInPrediction=FALSE
)

swissRes$summary

plot(swissRes$predict[["predict"]], main="predicted rain") 
plot(swissBorder, add=TRUE)

}

 