gmrfPrecUncond = function(x, 
			N = attributes(x)$Nx, Ny=attributes(x)$Ny,
			params = attributes(x)$model,
			border=params["kappa"]+1){
		if(!length(N) | !length(Ny) | !length(border) | !length(params))	
			warning("border, N, Ny, or params were not supplied")
		if(border < 1) warning("border should be >1")

		
		if(!all(c("scale","prec","kappa","cellSize") %in% names(params))) {
			warning("params needs scale, prec, kappa, cellSize")
			print(params)
		}
		
	Nx=N
		
		
	Ncells = Nx*Ny
	
	topCoords = rep(1:Nx, border) + 1i * rep(1:border, rep(Nx, border))
	topCells = Re(topCoords) + (Im(topCoords)-1)*Nx
	
	bottomCoords = rep(1:Nx, border) + 1i * rep(seq(Ny-border+1, Ny), rep(Nx, border))
	bottomCells = Re(bottomCoords) + (Im(bottomCoords)-1)*Nx
	
	Nonside = Ny-2*border
	leftCoords = rep(1:border, Nonside) + 
			1i * rep(seq(border+1, len=Nonside), rep(border,Nonside))
	leftCells = Re(leftCoords) + (Im(leftCoords)-1)*Nx
	rightCoords = rep(seq(Nx-border+1, Nx), Nonside) + 
			1i * rep(seq(border+1, len=Nonside), rep(border,Nonside))
	rightCells = Re(rightCoords) + (Im(rightCoords)-1)*Nx
	
	allCells = c(topCoords, leftCoords, rightCoords, bottomCoords)
	distmat = abs(outer(allCells, allCells, FUN="-"))
	distmat[upper.tri(distmat)]=NA
#	distmat = dist(cbind(Re(allCells), Im(allCells)),diag=T)
	distmat = new("dsyMatrix", Dim=rep(length(allCells), 2),
			x=as.vector(distmat),#[lower.tri(distmat, diag=T)], 
			uplo="L")
	
	cellSize = params["cellSize"]
	scaleCell = params["scale"] * cellSize
	distmat = distmat*scaleCell
	covMat = distmat
	covMat@x = (2^(1-params["kappa"])/(params["prec"]*gamma(params["kappa"]))) * 
			covMat@x^params["kappa"] * besselK(covMat@x, nu=params["kappa"])
	diag(covMat) = 1/params["prec"]
#	covChol = chol(covMat)
#	covInvChol = solve(covChol)
#	precOuter = solve(covMat)
	
	allCells = c(topCells, leftCells, rightCells, bottomCells)
	InnerPrecision = x[-allCells, -allCells]
	InnerPrecInvChol = solve(chol(InnerPrecision))
	
	
	
	A = x[allCells,-allCells]
	Aic = A %*% InnerPrecInvChol
	AQinvA = Aic %*% t(Aic)
	

	precOuter = solve(covMat) + AQinvA

	x[allCells,allCells] = precOuter
	x
}
