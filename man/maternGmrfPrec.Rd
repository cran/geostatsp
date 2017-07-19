\name{maternGmrfPrec}
\alias{maternGmrfPrec}
\alias{maternGmrfPrec.default}
\alias{maternGmrfPrec.dgCMatrix}
\alias{NNmat}
\alias{NNmat.default}
\alias{NNmat.Raster}
\title{Precision matrix for a Matern spatial correlation}
\description{
Produces the precision matrix for a Gaussian random field on a regular square lattice, using a Markov random field approximation. 
}

\usage{ 
maternGmrfPrec(N, ...)
\method{maternGmrfPrec}{dgCMatrix}(N, 
	param=c(variance=1, range=1, shape=1, cellSize=1),
  adjustEdges=FALSE,...) 
\method{maternGmrfPrec}{default}(N, Ny=N, 	
  param=c(variance=1, range=1, shape=1, cellSize=1),
  adjustEdges=FALSE, ...)
NNmat(N, Ny=N, nearest=3, adjustEdges=FALSE)
\method{NNmat}{Raster}(N, Ny=N, nearest=3, adjustEdges=FALSE)
\method{NNmat}{default}(N, Ny=N, nearest=3, adjustEdges=FALSE)
}

\arguments{
\item{N}{Number of grid cells in the x direction, or a matrix denoting nearest neighbours.}
\item{Ny}{Grid cells in the y direction, defaults to \code{N} for a square grid}
\item{param}{Vector of model parameters, with named elements: \code{scale},
scale parameter for the correlation function; \code{prec}, precision parameter; \code{shape}, 
Matern differentiability parameter (0, 1, or 2); and \code{cellSize}, the size of the grid cells.  Optionally, 
\code{variance} and \code{range} can be given in place of \code{prec} and \code{scale}, when the former 
are present and the latter are missing the reciprocal of the former are taken.}
\item{adjustEdges}{If \code{TRUE}, adjust the precision matrix so it does not implicitly assume 
the field takes values of zero outside the specified region.  Defaults to \code{FALSE}.  Can be a character string specifying the parameters to use for the correction, such as \code{'optimal'} or \code{'optimalShape'}, with \code{TRUE} equivalent to \code{'theo'}}
\item{nearest}{Number of nearest neighbours to compute}
\item{...}{Additional arguments passed to \code{maternGmrfPrec.dsCMatrix}}
}
\value{
A sparse matrix \command{\link[Matrix]{dsCMatrix-class}} object, containing a precision matrix for a 
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
	matrix(NNmat(11, 11, nearest=5)[,11*5+6],11, 11)

	params=c(range = 3,	shape=2, variance=5^2)
	
	myGrid = squareRaster(extent(0,20,0,10), 40)
		
	# precision matrix without adjusting for edge effects
	precMat =maternGmrfPrec(N=myGrid, param=params) 
	
	attributes(precMat)$info$precisionEntries
	
	midcell = cellFromRowCol(myGrid, 
		round(nrow(myGrid)/2), round(ncol(myGrid)/2)) # the middle cell
	edgeCell = cellFromRowCol(myGrid, 5,5)# cell near corner

# show precision of middle cell 
	precMid=matrix(precMat[,midcell], 
		nrow(myGrid), ncol(myGrid), byrow=TRUE)

	precMid[round(nrow(precMid)/2) + seq(-5, 5), 
		round(ncol(precMid)/2) + seq(-3, 3)]

	# and with the adjustment
	precMatCorr =maternGmrfPrec(
		N = myGrid, param=params, 
		adjustEdges=TRUE) 

if(interactive() | Sys.info()['user'] =='patrick') {	

# variance matrices
	varMat = Matrix::solve(precMat)
	varMatCorr = Matrix::solve(precMatCorr)

# compare covariance matrix to the matern
	xseq = seq(-ymax(myGrid), ymax(myGrid), len=1000)/1.5
	plot(xseq, matern(xseq, param=params),
	 type = 'l',ylab='cov', xlab='dist',
	 ylim=c(0, params["variance"]*1.1),
	 main="matern v gmrf")

	# middle cell
	varMid=matrix(varMat[,midcell], 
		nrow(myGrid), ncol(myGrid), byrow=TRUE)
	varMidCorr=matrix(varMatCorr[,midcell], 
		nrow(myGrid), ncol(myGrid), byrow=TRUE)
	xseqMid = yFromRow(myGrid) - yFromCell(myGrid, midcell)	
	points(xseqMid, varMid[,colFromCell(myGrid, midcell)], 
		col='red')
	points(xseqMid, varMidCorr[,colFromCell(myGrid, midcell)],
		 col='blue', cex=0.5)

	# edge cells
	varEdge=matrix(varMat[,edgeCell], 
	  nrow(myGrid), ncol(myGrid), byrow=TRUE)
	varEdgeCorr = matrix(varMatCorr[,edgeCell], 
	  nrow(myGrid), ncol(myGrid), byrow=TRUE)
	xseqEdge = yFromRow(myGrid) - yFromCell(myGrid, edgeCell)
	points(xseqEdge, 
		varEdge[,colFromCell(myGrid, edgeCell)], 
		pch=3,col='red')
	points(xseqEdge, 
	  varEdgeCorr[,colFromCell(myGrid, edgeCell)], 
	  pch=3, col='blue')
	
	legend("topright", lty=c(1, NA, NA, NA, NA), 
	  pch=c(NA, 1, 3, 16, 16),
		col=c('black','black','black','red','blue'),
		legend=c('matern', 'middle','edge','unadj', 'adj')
		)


	# construct matern variance matrix

	myraster = attributes(precMat)$raster
	covMatMatern = matern(myraster, param=params)
 
 	prodUncor = crossprod(covMatMatern, precMat)
 	prodCor = crossprod(covMatMatern, precMatCorr)

 	quantile(Matrix::diag(prodUncor),na.rm=TRUE)
 	quantile(Matrix::diag(prodCor),na.rm=TRUE)
 	
 	quantile(prodUncor[lower.tri(prodUncor,diag=FALSE)],na.rm=TRUE)	
 	quantile(prodCor[lower.tri(prodCor,diag=FALSE)],na.rm=TRUE)	
}
 	


}
