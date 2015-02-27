\name{profLlgm}
\alias{profLlgm}
\alias{informationLgm}
\title{
Joint confidence regions
}
\description{
Calculates profile likelihoods and approximate joint confidence regions
 for covariance parameters in linear geostatistical models.
}
\usage{
profLlgm(fit, mc.cores = 1, ...)
informationLgm(fit,  ...)
}

\arguments{
  \item{fit}{
Output from the \code{\link{lgm}} function
}
  \item{mc.cores}{
Passed to \code{\link[parallel]{mcmapply}}
}
  \item{\dots}{
For \code{profLlgm}, one or more vectors of parameter values 
at which the profile likelihood will be calculated, with names 
corresponding to elements of \code{fit$param}.  For \code{informationLgm},
arguments passed to  \code{\link[numDeriv]{hessian}}
}
}

\value{
\item{one or more vectors}{of parameter values}
\item{logL}{
A vector, matrix, or multi-dimensional array of profile likelihood values for every 
 combination of parameter values supplied.
}
\item{full}{Data frame with profile likelihood values and estimates of 
model parameters}
\item{prob,breaks}{vector of probabilities and chi-squared derived 
likelihood values associated with those probabilities}
\item{MLE,maxLogL}{Maximum Likelihood Estimates of parameters and 
log likelihood evaluated at these values}
\item{basepars}{combination of starting values for parameters re-estimated for each
profile likelihood and values of parameters which are fixed.}
\item{col}{vector of colours with one element fewer than the number of probabilities}
\item{ci,ciLong}{when only one parameter is varying, a matrix of confidence intervals (in both wide and long format)
is returned.}
}

\author{
Patrick Brown
}



\seealso{
\code{\link{lgm}}, \code{\link[parallel]{mcmapply}}, \code{\link[numDeriv]{hessian}} 
}
\examples{

# this example is time consuming
# the following 'if' statement ensures the CRAN
# computer doesn't run it
if(interactive() | Sys.info()['user'] =='patrick') {

library('geostatsp')
data('swissRain')

swissFit = lgm(data=swissRain, formula=rain~ SRTM_1km,
		grid=10, covariates=swissAltitude,
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE,reml=TRUE,
		param=c(anisoAngleDegrees=37,anisoRatio=7.5,
		range=50000))


x=profLlgm(swissFit,
		anisoAngleDegrees=seq(30, 43 , len=4)
)


plot(x[[1]],x[[2]], xlab=names(x)[1],
		ylab='log L',
		ylim=c(min(x[[2]]),x$maxLogL),
		type='n')
abline(h=x$breaks[-1],
		col=x$col,
		lwd=1.5)
axis(2,at=x$breaks,labels=x$prob,line=-1.2,
	tick=FALSE,
		las=1,padj=1.2,hadj=0)
abline(v=x$ciLong$par,
		lty=2,
		col=x$col[as.character(x$ciLong$prob)])
lines(x[[1]],x[[2]])



}

}

