\name{matern}
\alias{matern}
\alias{matern.default}
\alias{matern.dist}
\alias{matern.Raster}
\alias{matern.SpatialPoints}
\alias{matern.SpatialPointsDataFrame}

\title{Evaluate the Matern correlation function}
\description{
Returns the Matern covariance for the distances supplied.
}
\usage{
	matern( x, y=NULL, param=c(range=1, variance=1, shape=1))
	\method{matern}{SpatialPoints}(x,  y=NULL, param)
	\method{matern}{default}( x, y=NULL, param)
	\method{matern}{dist}( x, y=NULL, param)
	\method{matern}{Raster}( x, y=NULL, param)
	\method{matern}{SpatialPointsDataFrame}(x, y=NULL, param)
}

\arguments{
	\item{x}{A vector or matrix of distances, or \code{Raster} or
	 \code{SpatialPoints} of locations, see Details below.}
	 \item{y}{Covariance is calculated for the distance between locations in 
	 \code{x} and \code{y}. If \code{y=NULL}, covariance 
	 of \code{x} with itself is produced. 
	 However, if \code{x} is a matrix or vector it is assumed to 
	 be a set of distances and \code{y} is ignored.}
	\item{param}{A vector of named model parameters with, at a minimum names
 \code{range} and  \code{shape} (see Details), and optionally \code{variance} (defaults to 1).  
	For Geometric Anisotropy add 
	\code{anisoRatio} and either \code{anisoAngleDegrees} or  
		\code{anisoAngleRadians}
	}
}

\value{
	When \code{x} is a vector or matrix or object of class \code{dist}, a vector or matrix
	of covariances is returned.  
	With \code{x} being \code{SpatialPoints},  \code{y} must also be \code{SpatialPoints} and 
	a matrix of correlations between \code{x} and \code{y}
	is returned.	
	When \code{x} is a Raster, and \code{y} is a single location 
	a Raster of covariances between each pixel centre of \code{x} and \code{y} is returned.
}

\details{		
The formula for the Matern correlation function is
\deqn{M(x) =  \frac{variance}{\Gamma(shape)} 2^{shape-1} 
x (\sqrt{8 shape}/ range)^{shape}  
besselK(x \sqrt{8 shape}/ range, shape) }{%
             M(x) =   (variance / Gamma(shape)) 2^(shape-1) 
x [sqrt(8 shape)/ range]^shape  
besselK[x sqrt(8 shape)/ range, shape] 
}

The \code{range} argument is sqrt(8*shape)*phi.geoR, sqrt(8*shape)*scale.whittle.RandomFields, and
2*scale.matern.RandomFields.

Geometric anisotropy is only available when 
\code{x} is a \code{Raster} or \code{SpatialPoints}.  The parameter 'anisoAngle' refers to 
rotation of the coordinates anti-clockwise by the specified amount prior to 
calculating distances, which has the effect that the contours of the correlation function
are rotated clockwise by this amount. \code{anisoRatio} is the amount the Y coordinates are 
divided by
by post rotation prior to calculating distances.  A large value of \code{anisoRatio}
makes the Y coordinates smaller and increases the correlation in the
Y direction.

When \code{x} or \code{y} are rasters, cells are indexed row-wise 
starting at the top left.

}

\examples{
# example with raster
myraster = raster(nrows=40,ncols=60,xmn=-3,xmx=3,ymn=-2,ymx=2)
param = c(range=2, shape=2,	anisoRatio=2, anisoAngleDegrees=-25)

# plot correlation of each cell with the origin
myMatern = matern(myraster, c(0,0), param=param)


plot(myMatern, main="anisortopic matern")


# correlation matrix for all cells with each other
myraster = raster(nrows=4,ncols=6,xmn=-3,xmx=3,ymn=-2,ymx=2)
myMatern = matern(myraster, param=c(range=2, shape=2))
dim(myMatern)

# plot the cell ID's
values(myraster) = seq(1, ncell(myraster))
mydf = as.data.frame(myraster, xy=TRUE)
plot(mydf$x, mydf$y, type='n', main="cell ID's")
text(mydf$x, mydf$y, mydf$layer)
# correlation between bottom-right cell and top right cell is
myMatern[6,24]

# example with points
mypoints = SpatialPointsDataFrame(cbind(runif(5), runif(5)),data=data.frame(id=1:5))
matern(mypoints, param=param)

 
# example with vector of distances
range=3
distVec = seq(0, 2*range, len=100)
shapeSeq = c(0.5, 1, 2,20)
theCov = NULL
for(D in shapeSeq) {
	theCov = cbind(theCov, matern(distVec, param=c(range=range, shape=D)))
}
matplot(distVec, theCov, type='l', lty=1, xlab='distance', ylab='correlation',
	main="matern correlations")
legend("topright", fill=1:length(shapeSeq), legend=shapeSeq,title='shape')
# exponential

distVec2 = seq(0, max(distVec), len=20)
points(distVec2, exp(-2*(distVec2/range)),cex=0.5)
# gaussian
points(distVec2, exp(-2*(distVec2/range)^2), col='blue',cex=0.5)
legend("right", pch=1, col=c('black','blue'), legend=c('exp','gau'))

# comparing to geoR and RandomFields

covGeoR = covRandomFields = NULL

for(D in shapeSeq) {
	covGeoR = cbind(covGeoR, 
		geoR::matern(distVec, phi=range/sqrt(8*D), kappa=D))
	covRandomFields = cbind(covRandomFields,
		RandomFields::CovarianceFct(distVec, model="matern", 
			param = c(mean=0, variance=1, nugget=0, 
			scale=range/2, nu=D) ))
}

matplot(distVec, theCov, type='l', lty=1, xlab='distance', ylab='correlation',
	main="geoR, RandF, and geostatsp")
matpoints(distVec, covGeoR, cex=0.5, pch=1)
matpoints(distVec, covRandomFields, cex=0.5, pch=2)

legend("topright", lty=c(1,NA,NA), pch=c(NA, 1, 2), 
	legend=c("geostatsp","geoR","RandomFields"))

 
}


