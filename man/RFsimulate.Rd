\name{RFsimulate}
\alias{RFsimulate}
\alias{modelRandomFields}

\title{Simulation of Random Fields}

\description{
  This function simulates conditional and unconditional Gaussian random fields:
  
  Here, only the simulation of Gaussian random fields is described.
  For other kind of random fields (binary, max-stable, etc.) or
  more sophisticated approaches see \link[RandomFields]{RFsimulateAdvanced}.
}

\usage{
RFsimulate(model, x, y=NULL, z=NULL, T=NULL, grid,
 data, distances, dim, err.model, n=1, ...)
modelRandomFields(param, includeNugget=FALSE)
}

\arguments{
 \item{model}{object of class \code{\link[RandomFields]{RMmodel}}, a vector of named model
 parameters, or a matrix where each column is a model parameter
 }
 \item{x}{vector of x coordinates, or object
   which can be coerced to \code{\link[sp:GridTopology-class]{GridTopology}}.}
 \item{y}{optional vector of y coordinates}
 \item{z}{optional vector of z coordinates}
 \item{T}{optional vector of time coordinates, 
 \code{T} must always be an equidistant vector.
 Instead of \code{T=seq(from=From, by=By, len=Len)} one may also write
 \code{T=c(From, By, Len)}.
 }
 \item{grid}{logical; \code{RandomFields} can
 find itself the correct value in nearly all cases.
 See also \link[RandomFields]{RFsimulateAdvanced}.
 }
 \item{data}{For conditional simulation and random imputing only.
 If \code{data} is missing, unconditional
 simulation is performed.Object of class
 \command{\link[sp]{SpatialPointsDataFrame}};
 coordinates and response values of
 measurements in case that conditional simulation is to
 be performed
 }
 \item{distances}{another alternative to pass the (relative)
 coordinates.
 }
 \item{dim}{
 Only used if \code{distances} are given.
 }
 \item{err.model}{For conditional simulation and random imputing only.\cr
 Usually \code{err.model=RMnugget(var=var)}, or not given at all
 (error-free measurements).
 }
 \item{n}{number of realizations to generate.
 }
 \item{...}{for advanced use:
 further options and control parameters for the simulation
 that are passed to and processed by \command{\link[RandomFields]{RFoptions}}
 }
 \item{param}{A vector of named parameters}
 \item{includeNugget}{If \code{FALSE}, the \code{nugget} parameter is ignored. }
}

\details{
If \code{model} is a matrix, a different set of parameters is used for each simulation.  If 
\code{data} has the same number of columns as \code{model} has rows, 
a different column \code{i} is used with parameters in row \code{i}.   
}

\value{
If \code{x} is a raster, \code{RFsimulate} returns a \code{RasterLayer} or 
\code{RasterBrick} is returned.  Otherwise a
\code{SpatialPointsDataFrame} or \code{SpatialGridDataFrame} .  \code{modelRandomFields} returns
an object of class   \command{\link[RandomFields]{RMmodel}}.
 }


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}, Patrick Brown \email{patrick.brown@utoronto.ca}
 \url{http://pbrown.ca}
}

\seealso{
 \command{\link[RandomFields]{RFsimulate}},
 \command{\link[RandomFields]{RFempiricalvariogram}},
 \command{\link[RandomFields]{RFfit}},
 \command{\link[RandomFields]{RFgetModelInfo}},
 \command{\link[RandomFields]{RFgui}},
 \command{\link[RandomFields]{RMmodel}},
 \command{\link[RandomFields]{RFoptions}},
 \command{\link[RandomFields]{RFsimulateAdvanced}},
 \command{\link[RandomFields]{RFsimulate.more.examples}}
}

\examples{

model <- c(var=5, range=1,shape=0.5)

\dontrun{
	RandomFields::plot(modelRandomFields(model))
}

myraster = raster(nrows=60,ncols=60,xmn=0,ymn=0,xmx=10,ymx=10, 
		crs="+init=epsg:2081")

set.seed(0) 

simu <- geostatsp::RFsimulate(model, x=myraster, n=3)


plot(simu[[2]])



\dontrun{
data("swissRain")
swissRain$sqrtrain = sqrt(swissRain$rain)


# estimate parameters
swissRes =  lgm(data=swissRain, newdata=20, formula="sqrtrain",
	covariates=swissAltitude,  
	shape=1, fixShape=TRUE,
	aniso=FALSE, nugget=0, fixNugget=TRUE,
	nuggetInPrediction=FALSE
	)


# simulate from the random effect conditional on
#   the observed data

swissSim = geostatsp::RFsimulate(model=swissRes$param,
	data=swissRes[,'resid'],
	x=swissRes$predict,
	n=3
)

# plot the simulated random effect
plot(swissSim[[1]])
plot(swissBorder, add=TRUE)
}

}

\keyword{spatial}

