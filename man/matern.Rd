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
	matern( x, param=c(range=1, variance=1, shape=1),
	type=c('variance','cholesky','precision', 'inverseCholesky'),
	y=NULL)
	\method{matern}{SpatialPoints}(x,  param,
	type=c('variance','cholesky','precision', 'inverseCholesky'),
	y=NULL)
	\method{matern}{default}( x, param,
	type=c('variance','cholesky','precision', 'inverseCholesky'),
	y=NULL)
	\method{matern}{dist}( x, param,
	type=c('variance','cholesky','precision', 'inverseCholesky'),
	y=NULL)
	\method{matern}{Raster}( x,  param,
	type=c('variance','cholesky','precision', 'inverseCholesky'),
	y=NULL)
	\method{matern}{SpatialPointsDataFrame}(x, param,
	type=c('variance','cholesky','precision', 'inverseCholesky'),
	y=NULL)
}

\arguments{
	\item{x}{A vector or matrix of distances, or \code{Raster} or
	 \code{SpatialPoints} of locations, see Details below.}
	\item{param}{A vector of named model parameters with, at a minimum names
 \code{range} and  \code{shape} (see Details), and optionally \code{variance} (defaults to 1) and
 \code{nugget} (defaults to zero).  
	For Geometric Anisotropy add 
	\code{anisoRatio} and either \code{anisoAngleDegrees} or  
		\code{anisoAngleRadians}
	}
	\item{type}{specifies if the variance matrix, 
	the Cholesky decomposition of the variance matrix, 
	the precision matrix, or the inverse
	of the Cholesky L matrix is returned.}
\item{y}{Covariance is calculated for the distance between locations in 
	 \code{x} and \code{y}. If \code{y=NULL}, covariance 
	 of \code{x} with itself is produced. 
	 However, if \code{x} is a matrix or vector it is assumed to 
	 be a set of distances and \code{y} is ignored.}
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
 [x sqrt(8 shape)/ range]^shape  
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
param=c(shape=1,range=1,variance=1)
u=seq(0,4,len=20)
	uscale = sqrt(8*param['shape'])* u / param['range']
theMaterns = cbind(
	dist=u, 
	manual=	param['variance']*
			( 1/(gamma(param['shape']) * 2^(param['shape']-1)  ) ) * 
			uscale^param['shape'] * besselK( uscale , param['shape']),
	geostatsp=matern(u, param=param)
)
head(theMaterns)
matplot(theMaterns[,'dist'], 
	theMaterns[,c('manual','geostatsp')],
	col=c('red','blue'), type='l')
legend('topright', fill=c('red','blue'),
	legend=c('manual','geostatsp'))



# example with raster
myraster = raster(nrows=40,ncols=60,xmn=-3,xmx=3,ymn=-2,ymx=2)
param = c(range=2, shape=2,	anisoRatio=2, 
	anisoAngleDegrees=-25,variance=20)

# plot correlation of each cell with the origin
myMatern = matern(myraster, y=c(0,0), param=param)


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
mypoints = SpatialPoints(
	cbind(runif(8), runif(8))
	)
# variance matrix from points
m1=matern(mypoints, param=c(range=2,shape=1.4,variance=4,nugget=1))
# cholesky of variance from distances
c2=matern(dist(mypoints@coords), param=c(range=2,shape=1.4,variance=4,nugget=1),type='cholesky')

# check it's correct
quantile(m1- tcrossprod(c2))

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
legend("right", fill=1:length(shapeSeq), legend=shapeSeq,title='shape')
# exponential

distVec2 = seq(0, max(distVec), len=20)
points(distVec2, exp(-2*(distVec2/range)),cex=1.5, pch=5)
# gaussian
points(distVec2, exp(-2*(distVec2/range)^2), col='blue',cex=1.5, pch=5)
legend("bottomleft", pch=5, col=c('black','blue'), legend=c('exp','gau'))

# comparing to geoR and RandomFields

if (requireNamespace("RandomFields", quietly = TRUE) &
requireNamespace("geoR", quietly = TRUE) 
) { 

covGeoR = covRandomFields = NULL

for(D in shapeSeq) {
	covGeoR = cbind(covGeoR, 
		geoR::matern(distVec, phi=range/sqrt(8*D), kappa=D))
	covRandomFields = cbind(covRandomFields,
		RandomFields::RFcov(x=distVec, 
		model=RandomFields::RMmatern(nu=D, var=1,
				scale=range/2) ))
}



matpoints(distVec, covGeoR, cex=0.5, pch=1)
matpoints(distVec, covRandomFields, cex=0.5, pch=2)

legend("topright", lty=c(1,NA,NA), pch=c(NA, 1, 2), 
	legend=c("geostatsp","geoR","RandomFields"))
}
 
}


