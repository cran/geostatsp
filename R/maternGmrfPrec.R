getPrec = function(kappa,a)	{
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
		stop("shape parameter must be 0, 1 or 2")			
	}
	precEntries
}

maternGmrfPrec = function(N,...) {
	UseMethod("maternGmrfPrec")	
}

maternGmrfPrec.default = function(N, Ny=N,
		param=c(variance=1, range=1, shape=1, cellSize=1),
		 ...) {

	if(!'shape' %in% names(param)){
		warning("shape parameter appears to be missing")
	}
	
	theNNmat  = NNmat(N, Ny, nearest=param['shape']+1)

	maternGmrfPrec(theNNmat, param,...)
	
}


maternGmrfPrec.matrix = function(N, ...) {
	
	N = as(N, "dsCMatrix")
	maternGmrfPrec(N, ...)
	
}

maternGmrfPrec.dsCMatrix = function(N, 	
		param=c(variance=1, range=1, shape=1, cellSize=1),
		adjustEdges=FALSE,adjustParam=FALSE,
		adjustShape=FALSE,
		adjustMarginalVariance=FALSE,...) {

	names(param) = gsub("^var$", "variance", names(param))
	if(!any(names(param)=='variance'))
		param['variance']=1
	if(!any(names(param)=='cellSize'))
		param['cellSize']=1
	
	paramInfo = list(
		original=param,
		theo = NULL,
		sameShape = NULL,
		sameRange=NULL,
		optimal = NULL,
		oneminusar = NULL
	)
	
	objfun = function(oparam,distVec,sqrtVar,range=NULL,shape=NULL){
		
		if(!is.null(range))
			oparam['range'] = as.vector(range)
		if(!is.null(shape))
			oparam['shape'] = as.vector(shape)
		
		oparam['variance']=1
		
		theM = geostatsp::matern(
				x=distVec, 
				param=oparam,
		)
		sum( (sqrt(theM) - sqrtVar)^2)
	}


#	data('nn32')
#	tryprec = nn64
	tryprec = N
	Nx = attributes(tryprec)$Nx
	Ny = attributes(tryprec)$Ny
	
	midcellCoord = c(round(Nx*0.4),round(Ny*0.4)) # the middle cell
	midcell = c(Nx*(Ny-midcellCoord[2]) + midcellCoord[1])
	midVec = sparseMatrix(midcell,1,x=1,
			dims=c(ncol(N),1))

	distVecFull =  expand.grid(x=seq(1, Nx),
			y=seq(Ny,1))
	buffer = param['shape']+2
	
	cellVec = seq(1, nrow(distVecFull)) 
	
	isInner =  
			distVecFull[,'x']> buffer & 
					distVecFull[,'y']>buffer &
							distVecFull[,'x']< (Nx-buffer) &
							distVecFull[,'y'] < (Ny-buffer) 

	distVecFull = t(distVecFull)-midcellCoord 
	distVecFull = sqrt(apply(distVecFull^2,2,sum))	
		
	
	
	# if 1 - AR parameter supplied
	if(all(c("oneminusar","shape") %in% names(param))){
		a = 4/(1-param['oneminusar'])	
		

		paramInfo$theo = c(param[c('shape','cellSize','variance')],
				rangeInCells=as.numeric(sqrt(8*param['shape'])/sqrt(a-4)),
				a=as.vector(a)
			)
			
			
		if(min(c(Nx,Ny)<3*paramInfo$theo['rangeInCells'])){
				warning("grid is ", Nx, " by ", Ny,
						"which may be too small for range ",
						paramInfo$theo['rangeInCells'], " cells, 1- Ar param",
						param['oneminusar'])
			}
			
		paramInfo$theo['range'] = as.numeric(
				paramInfo$theo['rangeInCells']*param['cellSize']
			)
		paramInfo$oneminusar = param['oneminusar']
			
		precEntries=getPrec(param['shape'],a)
 
#	tryprec = nn64
		tryprec = N
		theN = tryprec@x
		theN = precEntries[theN]
		tryprec@x = theN
		
		varMid = solve(tryprec,midVec)
		varHere = varMid[midcell]
		
		whichDist = which(
				isInner &
				(
							distVecFull < (1.75*paramInfo$theo['rangeInCells'])
							) & 
						( distVecFull > 0))
		
		distVec = distVecFull[whichDist]
		varMid = varMid[whichDist]
		
		
		paramInfo$empirical = data.frame(x=c(0,distVec),
				y=c(varHere,varMid))			
		
		varMid = varMid/varHere
		
		sqrtVar = sqrt(varMid)
				

		
		startparam = c(param['shape'],
				range=as.vector(paramInfo$theo['rangeInCells'])
		)
		newPar = optim(
				startparam, objfun,
				lower=startparam/4,
				upper=startparam*4,
				distVec=distVec,sqrtVar=sqrtVar,
				method='L-BFGS-B'
		)
		paramInfo$optimal = c(
				newPar$par['shape'],
				param[c('variance','cellSize')],
				rangeInCells=as.vector(newPar$par['range']),
				newPar$par['range']*param['cellSize']
		)

		startparam=paramInfo$optimal['rangeInCells']
		names(startparam)='range'
		optSameShape = optim(
				startparam, objfun,
				lower=startparam/4,
				upper=startparam*4,
				shape=param['shape'],
				distVec=distVec,sqrtVar=sqrtVar,
				method='L-BFGS-B'
		)$par
		
		paramInfo$sameShape = paramInfo$theo
		paramInfo$sameShape['rangeInCells']=
				as.vector(optSameShape)
		paramInfo$sameShape['range']=
				as.vector(optSameShape*param['cellSize'])


		
		
	} else if(all(c('range','shape') %in% names(param))){
		param['rangeInCells'] = as.numeric(param['range']/param['cellSize'])

		scale = as.numeric(
				sqrt(8*param['shape'])/param["rangeInCells"]
		)
		
		
		a = (scale^2 + 4) 

		precEntries=getPrec(param['shape'],a)
#	data('nn32')
#	tryprec = nn64
		tryprec = N
		theN = tryprec@x
		theN = precEntries[theN]
		tryprec@x = theN
		
		varMid = solve(tryprec,midVec)
		varHere = varMid[midcell]
		
		
		
 
		
		
		if(min(c(Nx,Ny)<3*param['rangeInCells'])){
			warning("grid is ", Nx, " by ", Ny,
					"which may be too small for range ",
					paramInfo$theo['rangeInCells'], " cells")
		}
		

		
		whichDist = which(isInner &
				(
							distVecFull < (1.75*param['rangeInCells'])
							) & 
						( distVecFull > 0))
		
		distVec = distVecFull[whichDist]
		varMid = varMid[whichDist]/varHere
		sqrtVar = sqrt(varMid)

		
		paramInfo$theo = c(param, 	a=as.vector(a))
		
		if(adjustParam){
			startparam = c(#param['shape'],
				range=as.vector(param['rangeInCells'])
			)
			newPar = optim(
				startparam, objfun,
				lower=startparam/4,
				upper=startparam*4,
				distVec=distVec,sqrtVar=sqrtVar,
				shape=param['shape'],
				method='L-BFGS-B'
			)
		
			newRangeInCells = newPar$par['range']
			rangeCorrection = param['rangeInCells']/newRangeInCells
		
			newscale = sqrt(8*param['shape'])/(param['rangeInCells'] * 
					rangeCorrection)
			a = newscale^2 + 4
		
			precEntries = getPrec(param['shape'], a)	  

			tryprec = N
			theN = tryprec@x
			theN = precEntries[theN]
			tryprec@x = theN
		
			varMid = solve(tryprec,midVec)
			varHere = varMid[midcell]
					
			distVec = distVecFull[whichDist]
			
			varMid = varMid[whichDist]
			
					
			paramInfo$empirical = data.frame(x=c(0,distVec),
					y=c(varHere,varMid))			
			
			
			varMid = varMid/varHere
			
			sqrtVar = sqrt(varMid)

			paramInfo$theo['rangeInCells'] = as.vector(
					param['rangeInCells'] * rangeCorrection )
			paramInfo$theo['range'] = as.vector(
					paramInfo$theo['rangeInCells']*
							param['cellSize'])
			
		} 

		
		
		startparam = c(param['shape'],
				range=as.vector(param['rangeInCells']))

		newPar = optim(
				startparam, objfun,
				lower=startparam/4,
				upper=startparam*4,
				distVec=distVec,sqrtVar=sqrtVar,
				method='L-BFGS-B'
		)
		paramInfo$optimal = c(
				newPar$par['shape'],
				param[c('variance','cellSize')],
				rangeInCells=as.vector(newPar$par['range']),
				newPar$par['range']*param['cellSize']
		)
		
		startparam=paramInfo$optimal['rangeInCells']
		names(startparam)='range'
		optSameShape = optim(
				startparam, objfun,
				lower=startparam/4,
				upper=startparam*4,
				shape=param['shape'],
				distVec=distVec,sqrtVar=sqrtVar,
				method='L-BFGS-B'
		)$par
		
		paramInfo$sameShape = param 
		paramInfo$sameShape['rangeInCells']=
				as.vector(optSameShape)
		paramInfo$sameShape['range']=
				as.vector(optSameShape*param['cellSize'])
		
		
		startparam=paramInfo$optimal['shape']
		optSameRange = optim(
				startparam, objfun,
				lower=startparam/4,
				upper=startparam*4,
				range=param['rangeInCells'],
				distVec=distVec,sqrtVar=sqrtVar,
				method='L-BFGS-B'
		)$par
		paramInfo$sameRange = param 
		paramInfo$sameRange['shape']=
				as.vector(optSameRange)

		paramInfo$oneminusar = c(onemarinus=as.vector(1-4/a))
		
		
	} else {
		warning("param must have elements named shape  and either oneminusar or range")
		print(param)
	}

	if(adjustParam) {
		if(adjustShape){
			paramInfo$target = paramInfo$optimal
			} else { 
		paramInfo$target = paramInfo$sameShape
		}
	} else {
		paramInfo$target = paramInfo$theo
	}
	
	
	
# marginal precision
#	if(adjustParam){
if(adjustMarginalVariance){
	midQ = as(N[,midcell],'sparseVector')
		midQ@x = precEntries[midQ@x]
		

		paramHere = paramInfo$optimal
		
		paramHere = paramHere[
				c('shape','rangeInCells')]
		names(paramHere) = gsub("^rangeInCells$", "range",
				names(paramHere))
		paramHere['variance']=1
		midVar = matern(distVecFull[midQ@i],
				param= paramHere
						)	 
				
		marginalPrec =  sum(midQ@x * midVar)
		
		
} else {
		if(param['shape'] != 0) {
			  marginalPrec = (4*pi*param['shape'] *(a-4)^(param['shape'] ))
		  } else {
			  marginalPrec = (4*pi)
		  }
}

		# precEntries = precEntries/marginalPrec
 	#	precEntries = 
	#		precEntries/(marginalPrec*param['variance'])
#print(precEntries)
	precEntries = precEntries*exp(  
					 -log(marginalPrec) - 
					log(param['variance'])
	)
#	print(precEntries)
	paramInfo$empirical$yMult =
		paramInfo$empirical$y*exp(  
				log(marginalPrec) + 
						log(param['variance'])
	)
		
	
	theNNmat = N
	Nx=attributes(theNNmat)$Nx 
	Ny=attributes(theNNmat)$Ny 
	if(is.null(Nx)){
		Nx = Ny = sqrt(ncol(theNNmat))
		if(Nx!=round(Nx)){
			warning("N should have Nx and Ny attributes")
		}
	}
	theN = theNNmat@x
	theN = precEntries[theN]
	theNNmat@x = theN
		
		
	if(adjustEdges){
		
		distVecFull =  expand.grid(x=seq(1, Nx),
				y=seq(Ny,1))
		buffer = param['shape']+1
		
		cellVec = seq(1, nrow(distVecFull)) 
		
		whichInner = 
				distVecFull[,'x']> buffer & 
						distVecFull[,'y']>buffer &
						distVecFull[,'x']< (Nx-buffer) &
						distVecFull[,'y'] < (Ny-buffer) 
			
		innerCells = cellVec[whichInner]
		outerCells = cellVec[!whichInner]

		outerCoordsCartesian = SpatialPoints(
			distVecFull[outerCells,]*param['cellSize']
		)				
		if(adjustParam){
 				paramForM = paramInfo$target
		} else {
			paramForM = paramInfo$original	
			if(!any(names(paramForM)=='range')){
				paramForM= c(paramForM,
						paramInfo$theo['range']
						)
			}
		}
				
		covMat = matern(outerCoordsCartesian,
				param= paramForM)
			
		InnerPrecision = theNNmat[innerCells, innerCells]
			
			#A = x[allCells,-allCells]
			#InnerPrecInvChol = Matrix::solve(Matrix::chol(InnerPrecision))
			#Aic = A %*% InnerPrecInvChol
			# AQinvA = Aic %*% t(Aic)
		A = theNNmat[innerCells,outerCells]
		cholInnerPrec =Cholesky(InnerPrecision,LDL=FALSE,perm=TRUE)

		Aic = solve(cholInnerPrec, 
					solve(cholInnerPrec,A,system='P'),
					system='L')

		AQinvA = forceSymmetric(crossprod(Aic,Aic))

		covMatInv = Matrix::solve(covMat)
			
		precOuter = forceSymmetric(covMatInv + AQinvA)

		theNNmat[outerCells,outerCells] = precOuter
		theNNmat = forceSymmetric(theNNmat)
		paramInfo$adjust=c(edge=TRUE,
				param=adjustParam,
				shape=adjustShape)
		} else {
			paramInfo$adjust=c(edge=FALSE,
					param=adjustParam,
				shape=adjustShape)
		}
		paramInfo$adjust['marginalVariance'] = 
				adjustMarginalVariance
		paramInfo$precisionEntries = precEntries
		paramInfo$marginalPrec = as.vector(marginalPrec)
		
		theNNmat = drop0(theNNmat)

		
		attributes(theNNmat)$param=
				paramInfo

	attributes(theNNmat)$raster= list(
				nrows=Ny,
				ncols=Nx, 
				xmn=0,xmx=Nx*
						param['cellSize'],
				ymn=0, ymx=Ny*
						param['cellSize']
	)

					
	if(any(installed.packages()[,'Package'] == 'raster')) {
		attributes(theNNmat)$raster = 
				do.call(raster::raster,attributes(theNNmat)$raster)
	}
		
	return(theNNmat)
}
	
NNmat = function(N,Ny=N,nearest=3) {
		UseMethod("NNmat")	
}
	
NNmat.Raster = function(N,Ny=N,nearest=3) {
	NNmat(ncol(N),nrow(N),nearest)
}	

NNmat.default = function(N, Ny=N,nearest=3) {

	Nx = N
	Ncol = Nx
	Nrow=Ny
	Ncell = Nrow*Ncol
	
#	images.bresult = Matrix(data=0,nrow=Ncell, ncol=Ncell, sparse=T)
	Scell = 1:Ncell
	result = sparseMatrix(Scell, Scell, x=rep(1, Ncell),symmetric=TRUE)
#	diag(result) = 1
	

	
	# interior points
	oneN = c(1, -1, Ncol, -Ncol) # first neighbours up down right left
	if(nearest>=2) {
		twoN = c(Ncol-1, Ncol+1, -Ncol-1, -Ncol+1)  # first neighbours diagonal
		threeN = c(2, -2, 2*Ncol, -2*Ncol) # second neighbours up down right left
	} else {
		twoN = threeN = c()
	}
	
	fourN = c(3, -3, 3i, -3i) # square
	fiveseq = c(-1, -2, 1, 2) # diagonals
	fiveN = c(outer(fiveseq, fiveseq*1i, FUN="+"))
	fiveN = fiveN[abs(Re(fiveN) )!= abs(Im(fiveN))]
	if(nearest>=3) {
		fourNindex = Re(fourN) + Im(fourN)*Ncol
		fiveNindex = Re(fiveN) + Im(fiveN)*Ncol
	} else {
		fourNindex = fiveNindex = c()
	}
	
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
	if(nearest==1) {
		whichTwo = theNcEntries==2
		theNcVec = theNcVec[whichTwo]
		theNcEntries = theNcEntries[whichTwo]
	}	
	if(nearest==2) {
		whichTwo = theNcEntries<=4
		theNcVec = theNcVec[whichTwo]
		theNcEntries = theNcEntries[whichTwo]
	}	
	

			
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
	
	result = forceSymmetric(result)
	
	attributes(result)$Nx = Nx
	attributes(result)$Ny = Ny

	
	return(result)
}