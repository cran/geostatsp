\name{likfit}
\alias{likfit}
\alias{likfit.default}
\alias{likfit.SpatialPointsDataFrame}


\title{Likelihood Based Parameter Estimation for Gaussian Random Fields}

\description{
  \emph{Maximum likelihood} (ML) or \emph{restricted maximum likelihood} (REML)
  parameter estimation for (transformed) Gaussian random fields.
}



\usage{
likfit(geodata, ...)
\method{likfit}{SpatialPointsDataFrame}(geodata, formula, dist.rel.scale=100,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geodata}{
	An object of class \code{SpatialPointsDataFrame} or of a class suitable for 
	 \code{\link[geoR]{likfit}} in the \code{geoR} package.
}
\item{formula}{
A two sided formula, with the left side giving the response variable and the right side
the fixed effects variables
}
\item{dist.rel.scale}{the optimizer will work with a range parameter as a fraction of the diagonal distance of
the bounding box, with the default being 1/20 of the diagonal distance.  dist.rel.scale 
is the denominator of this fraction.
}
  \item{\dots}{
	additional arguments passed to 	 \code{\link[geoR]{likfit}} in the \code{geoR} package.
	
}
}

\value{
As \code{\link[geoR]{likfit}}  with additionally
  \item{beta.table}{a table of parameter estimates, standard errors, z scores, and p values}
}



\seealso{
\code{\link[geoR]{likfit}}
}
\examples{
n=100
mydat = SpatialPointsDataFrame(cbind(runif(n), runif(n)), 
	data=data.frame(cov1 = rnorm(n), cov2 = rpois(n, 0.5))
	)

# simulate a random field
mydat$U = GaussRF(mydat, model="whittle", param=c(mean=0, variance=1, nugget=0, 
				scale=0.4, alpha=2))

# add fixed effects
mydat$Y = -3 + 2*mydat$cov1 + 0.5*mydat$cov2 + rnorm(length(mydat), 0, 0.1)

myres = likfit(mydat, Y ~ cov1 + cov2, 
	cov.model="matern", kappa=1, fix.kappa=TRUE,
	ini.cov.pars = c(1, 0.2)
	)

myres$beta.table
}

