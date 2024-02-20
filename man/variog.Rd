\name{variog}
\alias{variog}
\alias{variog.default}
\alias{variog.SpatVector}
\alias{variogMcEnv}
\alias{variogMcEnv.default}
\alias{variogMcEnv.SpatVector}


\title{Compute Empirical Variograms and Permutation Envelopes}

\description{
These are wrappers for 
\Sexpr[results=rd]{c(
  '\\\\code{variog} in the \\\\code{geoR} package', 
  '\\\\command{\\\\link[geoR]{variog}}'
  )[1+requireNamespace('geoR', quietly=TRUE)]} 
 and 
 \Sexpr[results=rd]{c(
  '\\\\code{variog.mc.env} in the \\\\code{geoR} package', 
  '\\\\command{\\\\link[geoR]{variog.mc.env}}'
  )[1+requireNamespace('geoR', quietly=TRUE)]}.
}



\usage{
variog(geodata, \ldots)
\method{variog}{SpatVector}(geodata, formula, \ldots)
\method{variogMcEnv}{default}(geodata, \ldots)
\method{variogMcEnv}{SpatVector}(geodata, formula, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geodata}{
	An object of class \code{SpatVector} or of a class suitable for 
	\Sexpr[results=rd]{c(
  '\\\\code{variog} in the \\\\code{geoR} package', 
  '\\\\command{\\\\link[geoR]{variog}}'
  )[1+requireNamespace('geoR', quietly=TRUE)]} 
 }
	\item{formula}{A formula specifying the response variable and fixed effects portion of the model.
	The variogram is performed on the residuals.}
  \item{\dots}{
	additional arguments passed to 	\Sexpr[results=rd]{c(
  '\\\\code{variog} in the \\\\code{geoR} package', 
  '\\\\command{\\\\link[geoR]{variog}}'
  )[1+requireNamespace('geoR', quietly=TRUE)]} 
}
}

\value{
As \Sexpr[results=rd]{c(
  '\\\\code{variog} in the \\\\code{geoR} package', 
  '\\\\command{\\\\link[geoR]{variog}}'
  )[1+requireNamespace('geoR', quietly=TRUE)]} 
  and 
 \Sexpr[results=rd]{c(
  '\\\\code{variog.mc.env} in the \\\\code{geoR} package', 
  '\\\\command{\\\\link[geoR]{variog.mc.env}}'
  )[1+requireNamespace('geoR', quietly=TRUE)]} 
}



\seealso{
\Sexpr[results=rd]{c(
  '\\\\code{variog} in the \\\\code{geoR} package', 
  '\\\\command{\\\\link[geoR]{variog}}'
  )[1+requireNamespace('geoR', quietly=TRUE)] 
 } and 
 \Sexpr[results=rd]{c(
  '\\\\code{variog.mc.env} in the \\\\code{geoR} package', 
  '\\\\command{\\\\link[geoR]{variog.mc.env}}'
  )[1+requireNamespace('geoR', quietly=TRUE)]} 
}
\examples{
data("swissRain")
swissRain = unwrap(swissRain)
swissRain$lograin = log(swissRain$rain)
swissv= variog(swissRain, formula=lograin ~ 1,option="bin")
swissEnv = variogMcEnv(swissRain, lograin ~ 1, obj.var=swissv,nsim=9) 
if(!is.null(swissv)){
	plot(swissv, env=swissEnv, main = "Swiss variogram")
}
}


