library(geostatsp)
matrix(NNmat(7, 7)[,25], 7, 7)


params=c(range = 3,
		cellSize=0.5,
		shape=2,
		variance=5^2)

data("nn32")
# precision matrix without adjusting for edge effects
precMat =maternGmrfPrec(nn32, param=params) 
# and with the adjustment
precMatCorr =maternGmrfPrec(nn32, param=params, adjust.edges=TRUE) 

midcell = 32*16 + 16 # the middle cell
edgeCell = 32*5 + 5 # cell near corner

# show precision of cell 32,32 
precMid=matrix(precMat[,midcell], 32, 32, byrow=TRUE)
precMid[seq(16-4, 16+4), seq(16-4, 16+4)]

 
# variance matrices
	varMat = Matrix::solve(precMat)
	varMatCorr = Matrix::solve(precMatCorr)
	
	pdf("maternGmrfPred.pdf",height=5,width=7)
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
	dev.off()	
	
	# construct matern variance matrix
	Nx = attributes(precMat)$Nx
	Ny = attributes(precMat)$Ny
	cellSize = params["cellSize"]
	
	myraster = raster(nrows=Ny, ncols=Nx,
			xmn=0,ymn=0, xmx=Nx*cellSize, ymx=Ny*cellSize)
	covMatMatern = matern(myraster, param=params)
	
	prodUncor = covMatMatern %*% precMat
	prodCor = covMatMatern %*% precMatCorr
	
	
	
	quantile(Matrix::diag(prodUncor),na.rm=TRUE)
	quantile(Matrix::diag(prodCor),na.rm=TRUE)
	
	quantile(prodUncor[lower.tri(prodUncor,diag=FALSE)],na.rm=TRUE)	
	quantile(prodCor[lower.tri(prodCor,diag=FALSE)],na.rm=TRUE)	
	