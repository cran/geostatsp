
maternGmrfPrec = function(N,...) {
	UseMethod("maternGmrfPrec")	
}

maternGmrfPrec.default = function(N, Ny=N, ...) {


	theNNmat  = NNmat(N, Ny)

	maternGmrfPrec(theNNmat, ...)
	
}

maternGmrfPrec.matrix = function(N, ...) {
	
	N = as(N, "dgCMatrix")
	maternGmrfPrec(N, ...)
	
}

maternGmrfPrec.dgCMatrix = function(N, 	
		param=c(variance=1, range=1, rough=1, cellSize=1),
		adjust.edges=F,...) {

	names(param) = gsub("^var$", "variance", names(param))
	
	if(any(names(param)=="variance") & !any(names(param)=="prec"))
		param["prec"] = 1/param["variance"]
	if(any(names(param)=="range") & !any(names(param)=="scale"))
		param["scale"] = sqrt(8*param['rough'])/param["range"]

	if(!all( c("prec","scale","rough","cellSize")%in% names(param))) {
		warning("param must have elements named rough, cellSize, either prec or variance, and either scale or range")
	print(param)
	}
		
	theNNmat = N
	
	scale=param["scale"]
	prec=param["prec"]
	kappa=param["rough"]
	cellSize=param["cellSize"]
	
	scale = scale * cellSize
	a = (scale^2 + 4) 
	
 
	
	if(kappa == 0){
		precEntries = c(
				"1" = a,
				"2" = -1,
				"3" = 0,
				"4" = 0, 
				"5" =  0,
				"6" = 0)
	} else if(kappa==1) {	
		
		precEntries = c("1" = 4 + a^2,
					"2" = -2*a,
					"3" = 2,
					"4" = 1, 
					"5" =  0,
					"6" = 0)
		} else if(kappa==2) {
			precEntries = c("1" = a*(a*a+12),
					"2" = -3*(a*a+3),
					"3" = 6*a,
					"4" = 3*a, 
					"5" =  -3,
					"6" = -1)
		} else {
			stop("kappa must be 0, 1 or 2")			
		}
		
 
		
		if(kappa != 0) {
			marginalPrec = 4*pi*kappa*scale^(2*kappa)
		} else {
			marginalPrec = 4*pi
		}
 		precEntries = 
				precEntries*(prec /marginalPrec)
		
		theN = theNNmat@x
		theN = precEntries[theN]
		theNNmat@x = theN
		
	

		Nx=attributes(theNNmat)$Nx 
		Ny=attributes(theNNmat)$Ny 
		
		theNNmat = forceSymmetric(theNNmat)
		attributes(theNNmat)$model =param
		attributes(theNNmat)$Nx =Nx 
		attributes(theNNmat)$Ny =Ny
		
		if(adjust.edges){
			theNNmat = gmrfPrecUncond(theNNmat)
			attributes(theNNmat)$model =param
			attributes(theNNmat)$Nx =Nx 
			attributes(theNNmat)$Ny =Ny
		}
		
		
		return(theNNmat)
	}

NNmat = function(N, Ny=N) {

	Nx = N
	Ncol = Nx
	Nrow=Ny
	Ncell = Nrow*Ncol
	
#	images.bresult = Matrix(data=0,nrow=Ncell, ncol=Ncell, sparse=T)
	Scell = 1:Ncell
	result = sparseMatrix(Scell, Scell, x=rep(1, Ncell))
#	diag(result) = 1
	

	
	# interior points
	oneN = c(1, -1, Ncol, -Ncol) # first neighbours up down right left
	twoN = c(Ncol-1, Ncol+1, -Ncol-1, -Ncol+1)  # first neighbours diagonal
	threeN = c(2, -2, 2*Ncol, -2*Ncol) # second neighbours up down right left

	
	fourN = c(3, -3, 3i, -3i) # square
	fiveseq = c(-1, -2, 1, 2) # diagonals
	fiveN = c(outer(fiveseq, fiveseq*1i, FUN="+"))
	fiveN = fiveN[abs(Re(fiveN) )!= abs(Im(fiveN))]

	fourNindex = Re(fourN) + Im(fourN)*Ncol
	fiveNindex = Re(fiveN) + Im(fiveN)*Ncol
	
	if(all(c(Nrow,Ncol) >= 7)) {
		Scol = seq(4, Ncol-3)
		
		NeighbourIndexSeq = c(oneN, twoN, threeN, fourNindex,fiveNindex)
		NeighbourSeq = c(rep(2,length(oneN)), 
				rep(3 ,length(twoN)), 
				rep(4,length(threeN)), 
				rep(6 ,length(fourNindex)),
				rep(5 ,length(fiveNindex))
		)
		

	for(Drow in seq(4,Nrow-3)) {
		Prow =  (Drow-1)*Ncol
		for(Dcol in Scol){
 
			Dcell =  Prow + Dcol
			
			result[Dcell + NeighbourIndexSeq,Dcell] = NeighbourSeq
		}
	}
	}	

	# the borders
	
	theNc = cbind(oneN = c(1, -1, 1i, -1i), # first neighbours up down right left
	twoN = c(1-1i, -1+1i, -1i-1, 1i+1),  # first neighbours diagonal
	threeN = c(2, -2, 2i, -2i), # second neighbours up down right left
	fourN=fourN)	
	theNcVec = c(theNc[,"oneN"],
			theNc[,"twoN"],
			theNc[,"threeN"],
			theNc[,"fourN"],
			fiveN)
	theNcEntries = c(
			rep(2, dim(theNc)[1]),
			rep(3, dim(theNc)[1]),
			rep(4, dim(theNc)[1]),
			rep(6, dim(theNc)[1]),
			rep(5, length(fiveN))			
			)
	


			
	Scol = unique(c(1,2,3,Ncol-2,Ncol-1, Ncol))
	for(Drow in 1:Nrow) {
		Prow =  (Drow-1)*Ncol
		for(Dcol in Scol){
			Dcell = Prow + Dcol
			Ccell = Drow + 1i*Dcol
			Nhere = theNcVec + Ccell
			
			outsideBox = Re(Nhere) < 1 | Re(Nhere)>Nrow | 
					Im(Nhere) < 1 | Im(Nhere) > Ncol
			inBox = !outsideBox
			
			NhereIndex = (Re(Nhere) -1)*Ncol + Im(Nhere)
			
			result[NhereIndex[inBox],Dcell ] = theNcEntries[inBox]
		}			
	}
	Srow = unique(c(1,2,3,Nrow-2,Nrow-1, Nrow))
	Scol = seq(4, Ncol-3)
	for(Drow in Srow) {
		Prow =  (Drow-1)*Ncol
		for(Dcol in Scol){
			Dcell = Prow + Dcol
			Ccell = Drow + 1i*Dcol
			Nhere = theNcVec + Ccell
			
			outsideBox = Re(Nhere) < 1 | Re(Nhere)>Nrow | 
					Im(Nhere) < 1 | Im(Nhere) > Ncol
			inBox = !outsideBox
			
			NhereIndex = (Re(Nhere) -1)*Ncol + Im(Nhere)
			
			result[NhereIndex[inBox],Dcell ] = theNcEntries[inBox]
		}			
	}
	attributes(result)$Nx = Nx
	attributes(result)$Ny = Ny
	
	return(result)
}