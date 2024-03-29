\name{excProb}
\alias{excProb}
 

\title{
 Exceedance probabilities
}
\description{
Calculate exceedance probabilities pr(X > threshold) from a fitted geostatistical model.}
\usage{
excProb(x, threshold=0, random=FALSE, template=NULL, templateIdCol=NULL,
nuggetInPrediction=TRUE)
}
 
\arguments{
\item{x}{Output from either the \code{lgm} or \code{glgm} functions, or
a list of two-column matrices with columns named \code{x} and \code{y} containing the 
posterior distributions of random effects, as produced by \code{inla}.
}
\item{threshold}{ the value which the exceedance probability is calculated with respect to.
}
\item{random}{Calculate exceedances for the random effects, rather than the
predicted observations (including fixed effects). }
\item{template}{A \code{SpatRaster} or \code{SpatVector}
object which the results will be contained in.
}
\item{templateIdCol}{The data column in \code{template} corresponding to names of \code{marginals}
}
\item{nuggetInPrediction}{If \code{TRUE}, calculate exceedance probabilities of new observations by adding the 
nugget effect. Otherwise calculate probabilities for the latent process.  Ignored if \code{x} is output from \code{glgm}.
}
}
\details{
	When \code{x} is the output from \code{\link{lgm}}, pr(Y>threshold) is calculated using
	the Gaussian distribution using the Kriging mean and conditional variance.  When
	\code{x} is from the \code{glgm} function, 
	the marginal posteriors are numerically integrated to obtain pr(X > threshold).
}
\value{
Either a vector of exceedance probabilities or an object of the same class as \code{template}.
}

\examples{
	data('swissRain')
	swissRain = unwrap(swissRain)
	swissAltitude = unwrap(swissAltitude)
	swissBorder = unwrap(swissBorder)
	swissFit =  lgm("rain", swissRain, grid=30, 
		boxcox=0.5,fixBoxcox=TRUE,	covariates=swissAltitude)
	swissExc = excProb(swissFit, 20)
	mycol = c("green","yellow","orange","red")
	mybreaks = c(0, 0.2, 0.8, 0.9, 1)
	plot(swissBorder)
	plot(swissExc, breaks=mybreaks, col=mycol,add=TRUE,legend=FALSE)
	plot(swissBorder, add=TRUE)
	legend("topleft",legend=mybreaks, col=c(NA,mycol))

\donttest{

if(requireNamespace("INLA", quietly=TRUE) ) {
  INLA::inla.setOption(num.threads=2)
  # not all versions of INLA support blas.num.threads
  try(INLA::inla.setOption(blas.num.threads=2), silent=TRUE)

	swissRain$sqrtrain = sqrt(swissRain$rain)
	swissFit2 =  glgm(formula="sqrtrain",data=swissRain, grid=40, 
	covariates=swissAltitude,family="gaussian")
	swissExc = excProb(swissFit2, threshold=sqrt(30))
	swissExc = excProb(swissFit2$inla$marginals.random$space, 0,
		template=swissFit2$raster)
	
}

}

}




 	