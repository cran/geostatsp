\name{anisoMatrix}
\alias{anisoMatrix}

\title{
Calculates a 2-D matrix for geometric anisotropy
}
\description{
Calculates a 2-D matrix for geometric anisotropy, for use with \code{\link{GaussRF}}
}

\usage{
anisoMatrix(angle, range=NULL,scale=1/range)
}
\arguments{
  \item{angle}{
	angle of rotation, in radians
}
  \item{range}{
	range for x and y axes after rotation, ignored if scale is specified
}
  \item{scale}{
	scale for x and y axes after rotation
}
}
\details{
Produces matrix A such that
	(xnew, ynew)' = A (x,y)'	
}
\value{
A 2 by 2 matrix
}
\author{
Patrick Brown
}


\seealso{
\code{\link{GaussRF}}
}
\examples{

myAnisoMat = anisoMatrix(-pi/8, c(2, 0.5))

mymodel=list("$", var=2,  
		A=myAnisoMat,
		list("matern", nu=2))
		
myraster = raster(nrows=41,ncols=41,xmn=0,xmx=10,ymn=0,ymx=10)

myGRF=GaussRF(myraster, model=mymodel)

plot(myGRF)


}
