\name{variog}
\alias{variog}
\alias{variog.default}
\alias{variog.SpatialPointsDataFrame}
\alias{variog.mc.env}
\alias{variog.mc.env.default}
\alias{variog.mc.env.SpatialPointsDataFrame}


\title{Compute Empirical Variograms and Permutation Envelopes}

\description{
These are wrappers for \code{\link[geoR]{variog}} and \code{\link[geoR]{variog.mc.env}}
 in the \code{geoR} package.
}



\usage{
variog(geodata, \ldots)
\method{variog}{SpatialPointsDataFrame}(geodata, formula, \ldots)
\method{variog.mc.env}{default}(geodata, \ldots)
\method{variog.mc.env}{SpatialPointsDataFrame}(geodata, formula, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geodata}{
	An object of class \code{SpatialPointsDataFrame} or of a class suitable for 
	 \code{\link[geoR]{variog}} in the \code{geoR} package.}
	\item{formula}{A formula specifying the response variable and fixed effects portion of the model.
	The variogram is performed on the residuals.}
  \item{\dots}{
	additional arguments passed to 	 \code{\link[geoR]{variog}} in the \code{geoR} package.
}
}

\value{
As \code{\link[geoR]{variog}} or \code{\link[geoR]{variog.mc.env}} 
}



\seealso{
\code{\link[geoR]{variog}} and  \code{\link[geoR]{variog.mc.env}}.
}
\examples{
data(swissRain)
swissRain$lograin = log(swissRain$rain)
swissv= variog(swissRain, formula=lograin ~ 1,option="bin")
swissEnv = variog.mc.env(swissRain, lograin ~ 1, obj.var=swissv,nsim=9) 
plot(swissv, env=swissEnv)
}

