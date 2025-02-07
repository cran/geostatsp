\name{glgm-methods}
\docType{methods}
\alias{glgm}
\alias{glgm-methods}
\alias{glgm,ANY,ANY,ANY,ANY-method}
\alias{glgm,formula,SpatRaster,ANY,ANY-method}
\alias{glgm,formula,SpatVector,ANY,ANY-method}
\alias{glgm,formula,data.frame,SpatRaster,data.frame-method}


\alias{lgcp}

\title{
Generalized Linear Geostatistical Models
}
\description{
Fits a generalized linear geostatistical model or a log-Gaussian Cox process
using \code{inla} }
\usage{
\S4method{glgm}{ANY,ANY,ANY,ANY}(formula, data,  grid, covariates, buffer=0, shape=1, prior, ...) 
\S4method{glgm}{formula,SpatRaster,ANY,ANY}(formula, data,  grid, covariates, buffer=0, shape=1, prior, ...) 
\S4method{glgm}{formula,SpatVector,ANY,ANY}(formula, data,  grid, covariates, buffer=0, shape=1, prior, ...) 
\S4method{glgm}{formula,data.frame,SpatRaster,data.frame}(formula, data,  grid, covariates, buffer=0, shape=1, prior, ...) 
lgcp(formula=NULL, data,  grid, covariates=NULL, border, ...)
}
 
\arguments{
\item{data}{
An object of class \code{ SpatVector} containing the data.
}
\item{grid}{Either an integer giving the number of cells in the x direction, or a raster object which 
will be used for the spatial random effect.  If the cells in the raster are not square, the resolution in the y direction 
will be adjusted to make it so.}
\item{covariates}{Either a single raster, a list of rasters or a raster stack containing covariate values used when 
making spatial predictions.  Names of the raster layers or list elements correspond to names in the formula.  If
a covariate is missing from the data object it will be extracted from the rasters.  Defaults to \code{NULL} for an 
intercept-only model.
}
\item{formula}{Model formula, defaults to a linear combination of each of the layers in the \code{covariates} object.   
The spatial random effect should not be supplied but the default 
can be overridden with a 
 \code{ f(space,..)} term. For \code{glgm} the response variable defaults to the first variable in the \code{data} object, and
 \code{formula} can be an integer or character string specifying the response variable.  For \code{lgcp}, the formula
 should be one-sided.
}
\item{prior}{list with elements named \code{range}, \code{sd}, \code{sdObs}.  See Details.}
\item{shape}{Shape parameter for the Matern correlation function, must be 1 or 2.}
\item{buffer}{Extra space padded around the data bounding box to reduce edge effects.}
\item{border}{boundary of the region on which an LGCP is defined, passed to \code{\link[terra]{mask}}}
\item{...}{Additional options passed to  
%\Sexpr[results=rd]{c(
\code{inla} in the \code{INLA} package
%  '\\\\code{\\\\link[INLA]{inla}}'
%  )[1+requireNamespace('INLA', quietly=TRUE)]} 
}
}
\details{
This function performs Bayesian inference for generalized linear geostatistical models with INLA.  The Markov random field
approximation on a regular lattice is used for the spatial random effect.  The range parameter is the distance at which 
the correlation is 0.13, or
\deqn{cov[U(s+h), U(s)] = (2^{1-\nu}/Gamma(\nu))  d^\nu besselK(d, \nu) } 
\deqn{d= |h| \sqrt{8 \nu}/range}
where \eqn{\nu} is the shape parameter. The range parameter produced by \code{glgm} multiplies the range parameter from \code{INLA} by the cell size.

Elements of \code{prior} can be named \code{range}, \code{sd}, or \code{sdObs}.  Elements can consist of: 
\itemize{
\item a single value giving the prior median for penalized complexity priors (exponential on the sd or 1/range).
\item a vector \code{c(u=a, alpha=b)} giving an quantile and probability for pc priors.  For standard deviations alpha is an upper quantile, for the range parameter b = pr(1/range > 1/a).
\item a vector \code{c(lower=a, upper=b)} giving a 0.025 and 0.975 quantiles for the sd or range.
\item a list of the form \code{list(prior='loggamma', param=c(1,2))} passed directly to inla.
\item a two-column matrix of prior densities for the sd or range.
}


}

\value{
A list with two components named  \code{inla}, \code{raster}, and \code{parameters}.  \code{inla} contains the results of the call to the
\code{inla} function.  \code{raster} is a raster stack with the following layers:
\item{random.}{mean, sd, X0.0??quant: Posterior mean, standard deviation, and quantiles of the random effect}
\item{predict.}{mean, sd, X0.0??quant: same for linear predictors, on the link scale}
\item{predict.exp}{posterior mean of the exponential of the linear predictor}
\item{predict.invlogit}{Only supplied if a binomial response variable was used.}

\code{parameters} contains a list with elements:
\item{summary}{a table with parameter estimates and posterior quantiles}
\item{range, sd}{prior and posterior distributions of range and standard deviations}
}

\seealso{
%\Sexpr[results=rd]{c(
%  '\\\\code{inla} in the \\\\code{INLA} package', 
%  '\\\\code{\\\\link[INLA]{inla}}'
%  )[1+requireNamespace('INLA', quietly=TRUE)]}, 
  \url{https://www.r-inla.org}
}
\examples{
\dontrun{
# geostatistical model for the swiss rainfall data

if(requireNamespace("INLA", quietly=TRUE) ) {
  INLA::inla.setOption(num.threads=2)
  # not all versions of INLA support blas.num.threads
  try(INLA::inla.setOption(blas.num.threads=2), silent=TRUE)
} 

require("geostatsp")
data("swissRain")
swissRain = unwrap(swissRain)
swissAltitude = unwrap(swissAltitude)
swissBorder = unwrap(swissBorder)

swissRain$lograin = log(swissRain$rain)
swissFit =  glgm(formula="lograin", data=swissRain, 
	grid=30, 
	covariates=swissAltitude, family="gaussian", 
	buffer=2000,
	prior = list(sd=1, range=100*1000, sdObs = 2),
	control.inla = list(strategy='gaussian')
	)

if(!is.null(swissFit$parameters) ) {
	
	swissExc = excProb(swissFit, threshold=log(25))

	swissExcRE = excProb(swissFit$inla$marginals.random$space, 
		log(1.5),template=swissFit$raster)

	swissFit$parameters$summary

	matplot(
		swissFit$parameters$range$postK[,'x'],
		swissFit$parameters$range$postK[,c('y','prior')],
		type="l", lty=1, xlim = c(0, 1000),
		xlab = 'km', ylab='dens')
	legend('topright', lty=1, col=1:2, legend=c('post','prior'))

	plot(swissFit$raster[["predict.exp"]]) 

	mycol = c("green","yellow","orange","red")
	mybreaks = c(0, 0.2, 0.8, 0.95, 1)
	plot(swissBorder)
	plot(swissExc, breaks=mybreaks, col=mycol,add=TRUE,legend=FALSE)
	plot(swissBorder, add=TRUE)
	legend("topleft",legend=mybreaks, fill=c(NA,mycol))


	plot(swissBorder)
	plot(swissExcRE, breaks=mybreaks, col=mycol,add=TRUE,legend=FALSE)
	plot(swissBorder, add=TRUE)
	legend("topleft",legend=mybreaks, fill=c(NA,mycol))
}

		


# a log-Gaussian Cox process example

myPoints = vect(cbind(rbeta(100,2,2), rbeta(100,3,4)))


mycov = rast(matrix(rbinom(100, 1, 0.5), 10, 10), extent=ext(0, 1, 0, 1))
names(mycov)="x1"


if(requireNamespace("INLA", quietly=TRUE) ) {
  INLA::inla.setOption(num.threads=2)
  # not all versions of INLA support blas.num.threads
  try(INLA::inla.setOption(blas.num.threads=2), silent=TRUE)
}

res = lgcp(
	formula=~factor(x1),
	data=myPoints, 
	grid=squareRaster(ext(0,1,0,1), 20), covariates=mycov,
	prior=list(sd=c(0.9, 1.1), range=c(0.4, 0.41),
	control.inla = list(strategy='gaussian'), verbose=TRUE)
)
if(length(res$parameters)) {  
	plot(res$raster[["predict.exp"]])
	plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)
}
}
}


 	