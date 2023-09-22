\name{RFsimulate}
\alias{RFsimulate}
\alias{modelRandomFields}
\alias{RFsimulate}
\alias{RFsimulate-methods}
\alias{RFsimulate,ANY,SpatRaster-method}
\alias{RFsimulate,numeric,SpatRaster-method}
\alias{RFsimulate,numeric,SpatVector-method}
\alias{RFsimulate,RMmodel,SpatVector-method}
\alias{RFsimulate,RMmodel,SpatRaster-method}
\alias{RFsimulate,matrix,SpatRaster-method} 
\alias{RFsimulate,matrix,SpatVector-method} 
\alias{RFsimulate,data.frame,ANY-method} 


\title{Simulation of Random Fields}

\description{
  This function simulates conditional and unconditional Gaussian random fields, calling the function in the RandomFields package of the same name.
}

\usage{
\S4method{RFsimulate}{ANY,SpatRaster}(model, x,	data=NULL,
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{numeric,SpatRaster}(model, x,data=NULL,
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{numeric,SpatVector}(model, x, data=NULL, 
err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{RMmodel,SpatRaster}(model, x, data=NULL, 
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{RMmodel,SpatVector}(model, x, data=NULL, 
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{matrix,SpatRaster}(model, x, 	data=NULL, 
 err.model=NULL, n = nrow(model), ...)
\S4method{RFsimulate}{matrix,SpatVector}(model, x,	data=NULL, 
	err.model=NULL, n = nrow(model), ...)
\S4method{RFsimulate}{data.frame,ANY}(model, x,	data=NULL, 
	err.model=NULL, n = nrow(model), ...)
modelRandomFields(param, includeNugget=FALSE)
}

\arguments{
 \item{model}{object of class \code{\link[RandomFields]{RMmodel}}, a vector of named model
 parameters, or a matrix where each column is a model parameter
 }
 \item{x}{Object of type \code{\link[terra]{SpatRaster}} or \command{\link[terra]{SpatVector}}.}
 \item{data}{For conditional simulation and random imputing only.
 If \code{data} is missing, unconditional
 simulation is performed. Object of class
 \command{\link[terra]{SpatVector}};
 coordinates and response values of
 measurements in case that conditional simulation is to
 be performed
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
An object of the same class as \code{x}.
 }


\author{Patrick E. Brown \email{patrick.brown@utoronto.ca}
}

\seealso{
 \command{\link[RandomFields]{RFsimulate}},
 \command{\link[RandomFields]{RFfit}},
 \command{\link[RandomFields]{RFgetModelInfo}},
 \command{\link[RandomFields]{RFgui}},
 \command{\link[RandomFields]{RMmodel}},
 \command{\link[RandomFields]{RFoptions}},
 \command{\link[RandomFields]{RFsimulateAdvanced}},
 \command{\link[RandomFields]{RFsimulate.more.examples}}
}

\examples{
library('geostatsp')

# exclude this line to use the RandomFields package
options(useRandomFields = FALSE)

model1 <- c(var=5, range=1,shape=0.5)


myraster = rast(nrows=20,ncols=30,extent = ext(0,6,0,4), 
	crs="+proj=utm +zone=17 +datum=NAD27 +units=m +no_defs")

set.seed(0) 

simu <- RFsimulate(model1, x=myraster, n=3)

plot(simu[['sim2']])

xPoints = suppressWarnings(as.points(myraster))
# conditional simulation
firstSample = RFsimulate(
	c(model1, nugget=1), 
	x=xPoints[seq(1,ncell(myraster), len=100), ],
	n=3
)

secondSample = RFsimulate(
	model = cbind(var=5:3, range=seq(0.05, 0.25, len=3), shape=seq(0.5, 1.5, len=3)),
  err.model = 1,
	x= myraster,
	data=firstSample,n=4
)

plot(secondSample)



}

\keyword{spatial}

