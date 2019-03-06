\name{postExp}
\alias{postExp}


\title{Exponentiate posterior quantiles}
\description{
Converts a summary table for model parameters on the log scale to the natural or exponentiated scale.
}
\usage{
	postExp(x, 
		exclude = grep('^(range|aniso|shape|boxcox)', rownames(x)),
		invLogit=FALSE)
}

\arguments{
	\item{x}{a matrix or data frame as returned by \code{\link{glgm}} }
	\item{exclude}{vector of parameters not transformed, defaults to the range parameter}
	\item{invLogit}{Converts intercept parameter to inverse-logit scale when \code{TRUE}.  Can also be a vector of parameters to inverse-logit transform.}
}

\examples{
require("geostatsp")
data("swissRain")
swissRain$lograin = log(swissRain$rain)

if(requireNamespace('INLA', quietly=TRUE)) {
swissFit =  glgm(formula="lograin", data=swissRain, 
	grid=20, 
	covariates=swissAltitude/1000, family="gaussian", 
	prior = list(sd=1, range=100*1000, sdObs = 2),
	control.inla = list(strategy='gaussian', int.strategy='eb'),
	control.mode = list(theta=c(1.6542995, 0.7137123,2.2404179))
	)
postExp(swissFit$parameters$summary)

swissFitLgm =  lgm( formula="lograin",
	data=swissRain[1:60,], grid=20,
	covariates=swissAltitude/1000, 
	shape=1, fixShape=TRUE,
	aniso=FALSE, nugget=0, fixNugget=FALSE
)

swissFitLgm$summary
}

}
