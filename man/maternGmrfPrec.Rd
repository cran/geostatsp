\name{maternGmrfPrec}
\alias{maternGmrfPrec}
\alias{maternGmrfPrec.default}
\alias{maternGmrfPrec.matrix}
\alias{maternGmrfPrec.dsCMatrix}
\alias{NNmat}
\alias{NNmat.default}
\alias{NNmat.Raster}
\alias{gmrfPrecUncond}
\title{Precision matrix for a Matern spatial correlation}
\description{
Produces the precision matrix for a Gaussian random field on a regular square lattice, using a Markov random field approximation. 
}

\usage{ 
maternGmrfPrec(N, ...)
\method{maternGmrfPrec}{dsCMatrix}(N, 
	param=c(variance=1, range=1, shape=1, cellSize=1),
  adjustEdges=FALSE,...) 
\method{maternGmrfPrec}{matrix}(N, ...)
\method{maternGmrfPrec}{default}(N, Ny=N, 	
param=c(variance=1, range=1, shape=1, cellSize=1),
 ...)
NNmat(N, Ny=N, nearest=3)
\method{NNmat}{Raster}(N, Ny=N, nearest=3)
\method{NNmat}{default}(N, Ny=N, nearest=3)

gmrfPrecUncond(x, 
			N = attributes(x)$Nx, Ny=attributes(x)$Ny,
			param = attributes(x)$model,
			border=param["shape"]+1)
}

\arguments{
\item{N}{Number of grid cells in the x direction, or a matrix denoting nearest neighbours.}
\item{Ny}{Grid cells in the y direction, defaults to \code{N} for a square grid}
\item{param}{Vector of model parameters, with named elements: \code{scale},
scale parameter for the correlation function; \code{prec}, precision parameter; \code{shape}, 
Matern differentiability parameter (0, 1, or 2); and \code{cellSize}, the size of the grid cells.  Optionally, 
\code{variance} and \code{range} can be given in place of \code{prec} and \code{scale}, when the former 
are present and the latter are missing the reciprocal of the former are taken.}
\item{border}{number of cells from the edge of the grid to adjust for edge effects.}
\item{adjustEdges}{If \code{TRUE}, adjust the precision matrix so it does not implicitly assume 
the field takes values of zero outside the specified region.  Defaults to \code{FALSE}.  Can be a character string specifying the parameters to use for the correction, such as \code{'optimal'} or \code{'optimalShape'}, with \code{TRUE} equivalent to \code{'theo'}}

\item{x}{An MRF precision matrix produced by \code{maternGmrfPrec}}
\item{nearest}{Number of nearest neighbours to compute}
\item{...}{Additional arguments passed to \code{maternGmrfPrec.dgCMatrix}}
}
\value{
A sparse matrix \code{dgCMatrix} object, containing a precision matrix for a 
Gaussian random field or (from the \code{NNmat} function)  a matrix denoting neighbours.  

 }
\details{ The numbering of cells
is consistent with the \code{raster} package.  Cell 1 is the top left cell, with cell 2 being the cell to the right and numbering
continuing row-wise.

The 
nearest neighbour matrix \code{N} has: \code{N[i,j]=1} if \code{i=j}; 
takes a value 2 if \code{i} and \code{j} are first `rook' neighbours; 
3 if they are first `bishop' neighbours; 4 if they are second `rook' neighbours; 5 
if `knight' neighbours; and 6 if third `rook' neighbours.
\preformatted{
     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,]    0    0    0    6    0    0    0
[2,]    0    0    5    4    5    0    0
[3,]    0    5    3    2    3    5    0
[4,]    6    4    2    1    2    4    6
[5,]    0    5    3    2    3    5    0
[6,]    0    0    5    4    5    0    0
[7,]    0    0    0    6    0    0    0
}
}

\examples{
# produces the matrix above
	matrix(NNmat(7, 7)[,25], 7, 7)


	params=c(range = 3,
	cellSize=0.5,
	shape=2,
	variance=5^2)
	
	data("nn32")
	# precision matrix without adjusting for edge effects
	precMat =maternGmrfPrec(nn32, param=params) 
	
	if(interactive() | Sys.info()['user'] =='patrick') {	
	# and with the adjustment
	precMatCorr =maternGmrfPrec(nn32, param=params, adjustEdges=TRUE) 
	
	midcell = 32*16 + 16 # the middle cell
	edgeCell = 32*5 + 5 # cell near corner

# show precision of cell 32,32 
	precMid=matrix(precMat[,midcell], 32, 32, byrow=TRUE)
	precMid[seq(16-4, 16+4), seq(16-4, 16+4)]

# variance matrices
	varMat = Matrix::solve(precMat)
	varMatCorr = Matrix::solve(precMatCorr)

# compare covariance matrix to the matern
	xseq = seq(-20*params["cellSize"], 20*params["cellSize"], len=1000)
	plot(xseq, matern(xseq, param=params),
	 type = 'l',ylab='cov', xlab='dist',ylim=c(0, params["variance"]),
	 main="matern v gmrf")



	# middle cell
	varMid=matrix(varMat[,midcell], 32, 32, byrow=TRUE)
	varMidCorr=matrix(varMatCorr[,midcell], 32, 32, byrow=TRUE)
	xseqMid = params["cellSize"] *seq(-16,15) 	
	points(xseqMid, varMid[,16], col='red')
	points(xseqMid, varMidCorr[,16], col='blue', cex=0.5)

	# edge cells
	varEdge=matrix(varMat[,edgeCell], 32, 32, byrow=TRUE)
	varEdgeCorr=matrix(varMatCorr[,edgeCell], 32, 32, byrow=TRUE)
	xseqEdge = params["cellSize"] *seq(-5, 26)
	points(xseqEdge, varEdge[,5], pch=3,col='red')
	points(xseqEdge, varEdgeCorr[,5], pch=3, col='blue')
	
	legend("topright", lty=c(1, NA, NA, NA, NA), pch=c(NA, 1, 3, 16, 16),
		col=c('black','black','black','red','blue'),
		legend=c('matern', 'middle','edge','unadj', 'adj')
		)


	# construct matern variance matrix

	myraster = attributes(precMat)$raster
	covMatMatern = matern(myraster, param=params)
 
 	prodUncor = covMatMatern \%*\% precMat
 	prodCor = covMatMatern \%*\% precMatCorr



 	quantile(Matrix::diag(prodUncor),na.rm=TRUE)
 	quantile(Matrix::diag(prodCor),na.rm=TRUE)
 	
 	quantile(prodUncor[lower.tri(prodUncor,diag=FALSE)],na.rm=TRUE)	
 	quantile(prodCor[lower.tri(prodCor,diag=FALSE)],na.rm=TRUE)	
}
 	


}
