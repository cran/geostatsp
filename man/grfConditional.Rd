\name{grfConditional}
\alias{grfConditional}

\title{
Conditional simulation of a Gaussian Random Field
}
\description{
Simulates a Gaussian random field conditional on observed values.
}
\usage{
grfConditional(data,  model.fit, locations, Nsim, 
	covariates, rasterMethod = "ngb", fun, nugget.in.prediction=T) 
}
 
\arguments{
\item{data}{ A SpatialPointsDataFrame object containing the data to be conditioned upon.
}
  \item{model.fit}{
A fitted model produced by \code{\link{likfit}}
}
\item{locations}{A raster on which the field will be simulated.}
\item{Nsim}{Number of samples to simulate.}
  \item{covariates}{
The spatial covariates used in prediction, either a \code{\link[raster]{raster}} stack or list of rasters.
}
\item{rasterMethod}{Method for reprojecting or resampling rasters, either "bilinear" or "ngb".}
\item{fun}{A function applied to each realised surface and only the result saved.}
\item{nugget.in.prediction}{If TRUE, predict new observations by adding the nugget effect.  Otherwise
predict fitted values.  Only relevant for Box-Cox or log transformed data.}
}

\value{
\item{random}{the predictions and conditional variance of the random effects, on the same
raster as \code{locations}}
}


\seealso{
\code{\link[geoR]{krige.conv}}
}
\examples{
data(swissRain)


swissRain$elevation = extract(swissAltitude, swissRain)
swissRain$sqrtrain = sqrt(swissRain$rain)
swissFit =  likfit(swissRain, sqrtrain ~ elevation, 
	cov.model="matern", kappa = 2, fix.kappa=TRUE, 
	fix.psiA=FALSE, fix.psiR=FALSE,
	ini.cov.pars= c(0.6, 7000))

swissRaster = raster(extent(swissRain@bbox), ncols=20, nrows=20, 
	crs=swissRain@proj4string)	


swissSim = grfConditional(data=swissRain,  model.fit=swissFit, 
locations=swissRaster, Nsim=2, 
covariates=swissAltitude, nugget.in.prediction=FALSE) 

plot(swissSim[[2]]) 

# proportion the realisation above 200
swissSim2 = grfConditional(data=swissRain,  model.fit=swissFit, 
locations=swissRaster, Nsim=4, covariates=swissAltitude, 
fun = function(x) mean(x^2>200)  )

unlist(swissSim2) 
\dontrun{
load(url("http://www.filefactory.com/file/frd1mhownd9/n/CHE_adm0_RData"))
plot(spTransform(gadm, swissRain@proj4string), add=TRUE)
}

}

