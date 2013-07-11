\name{loglik.GRF}
\alias{loglik.GRF.SpatialPointsDataFrame}
\alias{loglik.GRF}
\alias{loglik.GRF.default}




\title{
Log-Likelihood for a Gaussian Random Field}
\description{
This function computes the value of the log-likelihood for a Gaussian random field.}
\usage{
loglik.GRF(geodata, ...)
% coords=geodata$coords, data=geodata$data,
%				trend = ~1, ini.cov.pars,
%				fix.nugget = FALSE, nugget = 0, 
%				fix.kappa = TRUE, kappa = 0.5, 
%				fix.lambda = TRUE, lambda = 1, 
%				fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, 
%				cov.model, realisations, lik.method = "ML",
%				components = TRUE, nospatial = TRUE, limits = pars.limits(), 
%				print.pars = FALSE, messages, ...)
\method{loglik.GRF}{SpatialPointsDataFrame}(geodata, formula, sigmasq, 
  phi, psiA = 1, psiR = 1, kappa = 1, lambda = 1, nugget = 0, ...)
}
\arguments{
  \item{geodata}{
A SpatialPointsDataFrame containing the data to be fit, or an object suitable for passing to
\code{\link[geoR]{loglik.GRF}} in \code{geoR}.  
}
  \item{formula}{
A formula object, specifies the response variable and covariates
}
  \item{sigmasq}{
Variance of the Gaussian random field
}
  \item{phi}{
Range of the GRF
}
  \item{psiA}{
Anisotropy rotation angle
}
  \item{psiR}{
Ratio of major to minor axis lengths for anisotropy, defaults to 1 for an isotropic process.
}
  \item{kappa}{
Roughness parameter for Matern correlation
}
  \item{lambda}{
Box-Cox parameter.
}
  \item{nugget}{
Observation-level variance
}
  \item{\dots}{
Additional arguments for \code{\link[geoR]{loglik.GRF}} in \code{geoR}. 
}
%\item{coords}{}
%\item{data}{}
%\item{	trend = ~1, ini.cov.pars,
%				fix.nugget = FALSE, nugget = 0, 
%				fix.kappa = TRUE, kappa = 0.5, 
%				fix.lambda = TRUE, lambda = 1, 
%				fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, 
%				cov.model, realisations, lik.method = "ML",
%				components = TRUE, nospatial = TRUE, limits = pars.limits(), 
%				print.pars = FALSE, messages, 
}
\details{
This function is essentially a wrapper for
\code{\link[geoR]{loglik.GRF}} in \code{geoR}.  
}
\value{
the value of the log likelihood.
}

\examples{
data(swissRain)
swissRain$lograin = log(swissRain$rain)
swissRain$elevation = extract(swissAltitude, swissRain)
loglik.GRF(swissRain, lograin ~ elevation, 
	cov.model="matern", kappa = 2, 
	phi=10000,sigmasq=0.5^2)

loglik.GRF(swissRain, lograin ~ elevation, 
	cov.model="matern", kappa = 2, 
	phi=5000,sigmasq=1)
}

