# needed for optimal parameters
if(FALSE) {
	
	theNN = NNmat(11, 11, nearest=5)
	theRaster = attributes(theNN)$raster
	matSub = matrix(theNN[,11*5+6],11, 11)
	cCell = c(6,6)
	
	numbn = table(matSub)
	numbn = numbn[setdiff(names(numbn), '0')]
	dput(numbn, '')
	
	ndist = c('1'=0, '2'=1, '3'=sqrt(2), '4'=2, '5'=sqrt(5), '6'=3, '7'=sqrt(10), '8'=sqrt(8), '9'=4) 
	
	dMat = spDists(
			xyFromCell(theRaster, 1:ncell(theRaster)),
			cbind(xFromCol(theRaster, cCell[1]), yFromRow(theRaster, cCell[2]))
	)
	values(theRaster) = dMat
	theIndex = raster(theRaster)
	values(theIndex) = as.vector(matSub)
	
	dVec = cbind(values(theIndex), values(theRaster))
	dVec = dVec[dVec[,1]> 0,]
	dVec = dVec[order(dVec[,1]),]  
	dVec = dVec[!duplicated(dVec[,1]),]
	rownames(dVec) = dVec[,1]
	ndist = dVec[,2]
	dput(ndist, '')
}

numbn = structure(c(1L, 4L, 4L, 4L, 8L, 4L, 4L, 8L, 4L, 8L, 8L, 4L), .Dim = 12L, .Dimnames = structure(list(
						matSub = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
								"11", "12")), .Names = "matSub"), class = "table")

ndist = structure(c(0, 1, 1.4142135623731, 2, 2.23606797749979, 3, 2.82842712474619, 
				3.16227766016838, 4, 3.60555127546399, 4.12310562561766, 5), .Names = c("1", 
				"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

nbcells = as.vector(outer((-4):4, 1i*(-4):4,"+"))
nbcells = nbcells[signif(Mod(nbcells),3)%in%signif(ndist,3)]
nbn = names(ndist)[match(signif(Mod(nbcells),3), signif(ndist,3))]
othercells = as.vector(outer(0:15, 1i*1:15,"+"))
othercells = othercells[Re(othercells) <= Im(othercells)]
otherdist = Mod(outer(othercells, nbcells, '-'))
dimnames(otherdist) = list(as.character(othercells), as.character(nbcells))


getPrec = function(shape,oneminusar)	{
	a = as.numeric((1-oneminusar)/4)
	logA = as.numeric(log(1-oneminusar) - log(4))
	# This is no longer lindgren et al's a!  it's 1/their a
	# and the conditional variance is different
	if(shape == 0){
		precEntries = c(
				"1" = 1,
				"2" = -a)
	} else if(shape==1) {			
		precEntries = c(
				"1" = 1+4*a^2,
				"2" = -2*a,
				"3" = 2*a^2,
				"4" = a^2, 
				"5" =  0,
				"6" = 0)
	} else if(shape==2) {
		precEntries = c(
				"1" = 1+12*a^2,
				"2" = -3*a-9*a^3,
				"3" = 6*a^2,
				"4" = 3*a^2, 
				"5" =  -3*a^3,
				"6" = -a^3)
		
	} else if(shape==3) {
		# shape 3 https://bitbucket.org/hrue/r-inla/src/ gmrflib/matern.c
#    val = SQR(SQR(a) + 6.0) + 12 * SQR(a);
		precEntries = c(
				"1" = (1+6*a^2)^2 + 12*a^2,
				"2" = -4*(a+9*a^3),
				"3" = 12*(a^2+2*a^4),
				"4" = 2*(3*a^2+8*a^4), 
				"5" = -12*a^3,
				"6" = -4*a^3,
				"7" = 6*a^4, #2,2
				"8" = 4*a^4, #3,1
				"9" = a^4
		)
	} else if(shape==4) {
		# values computed by yacas in gmrfShape2.R file
		precEntries = c(
				'1'=  180 * a^4 + 40 * a^2 + 1 ,
				'2'=  5 * (-20 * a^5 - 18 * a^3 - a) ,
				'3'=  20 * (6 * a^4 + a^2) ,
				'4'=  10 * (8 * a^4 + a^2) ,
				'5'=  10 * (-5 * a^5 - 3 * a^3) ,
				'6'=  5 * (-5 * a^5 - 2 * a^3) ,
				'7'=  30 * a^4 ,
				'8'=  20 * a^4 ,
				'9'=  5 * a^4 ,
				'10'=  -10 * a^5 ,
				'11'=  -5 * a^5 ,
				'12'=  -a^5  
		)
	} else if (shape == 5){
		precEntries = c(
#	  "1" = 400*a^6+540*a^4+60*a^2+1, 
				"1" = exp(log(400) + 6*logA)+exp(log(540)+4*logA)+exp(log(60)+2*logA)+1, 
				"2" = (-600)*a^5-180*a^3-6*a, 
				"3" = 30*(10*a^6+12*a^4+a^2), 
				"4" = 15*(15*a^6+16*a^4+a^2), 
				"5" = (-300)*a^5-60*a^3, 
				"6" = (-150)*a^5-20*a^3, 
				"7" = 30*(4*a^6+3*a^4), 
				"8" = 30*(3*a^6+2*a^4), 
				"9" = 3*(12*a^6+5*a^4), 
				"10" = (-60)*a^5, 
				"11" = (-30)*a^5, 
				"12" = (-6)*a^5, 
				"13" = 20*a^6, 
				"14" = 15*a^6, 
				"15" = 6*a^6, 
				"16" = a^6
		)
	} else {
		stop("shape parameter must be 0, 1, 2 or 3")			
	}
	precEntries
}

maternGmrfPrec = function(N,...) {
	UseMethod("maternGmrfPrec")	
}

maternGmrfPrec.default = function(N, Ny=N,
		param=c(variance=1, range=1, shape=1, cellSize=1),
		adjustEdges=FALSE,
		...) {
	
	if(!'shape' %in% names(param)){
		warning("shape parameter appears to be missing")
	}
	
	theNNmat  = NNmat(N, Ny, nearest=param['shape']+1, adjustEdges=adjustEdges)
	
	maternGmrfPrec(N=theNNmat, param=param, adjustEdges=adjustEdges, ...)
	
}


maternGmrfPrec.dgCMatrix = function(N,
		param=c(variance=1, range=1, shape=1, cellSize=1),
		adjustEdges=FALSE,
		...) {
	
	param['nugget']=0
	
	names(param) = gsub("^var$", "variance", names(param))
	if(!any(names(param)=='variance')) {		
		if ('conditionalVariance' %in% names(param)){
			param['variance']=NA
		} else {
			param['variance']=1
		}
	} else {
		if ('conditionalVariance' %in% names(param)){
			if(any(names(param)=='oneminusar')) {
				warning(
						"both conditionalVariance and variance were supplied, ignoring variance"
				)	
				param['variance']=NA
			} else {
				warning(
						"both conditionalVariance and variance were supplied, ignoring conditionalVariance"
				)	
				param['conditionalVariance']=NA
			}
		}
	}
	
	
	theraster = attributes(N)$raster
	if(is.null(theraster)) {
		if(!any(names(param)=='cellSize')){
			param['cellSize']=1
		}
		if(!all(c("Nx","Ny")%in% names(attributes(N)))) {
			Nx = Ny = sqrt(ncol(N))
			if(Nx!=round(Nx)){
				warning("N should have Nx and Ny attributes")
			} 
		} else {
			Nx = attributes(N)$Nx
			Ny = attributes(N)$Ny
		}
		theraster = list(
				nrows=Ny,
				ncols=Nx, 
				xmn=0,xmx=Nx*
						param['cellSize'],
				ymn=0, ymx=Ny*
						param['cellSize']
		)	
		
	} else { # the raster not null
		
		Nx = ncol(theraster)
		Ny = nrow(theraster)
		param['cellSize'] = 
				( attributes(theraster)$extent@xmax - 
					attributes(theraster)$extent@xmin) /
				Nx		
	}
	
	gmrfOrder = param['shape']+1
	
	adjustEdgesN = attributes(N)$adjustEdges
	if(is.null(adjustEdgesN)) adjustEdgesN = FALSE
	
	nearestN = attributes(N)$nearest
	if(is.null(nearestN)) nearestN = Inf
	
	adjustEdgesLogical = adjustEdges != FALSE
	
	# recompute N if the edge adjustment isn't appropriate
	if( (adjustEdgesLogical != adjustEdgesN) |
			(adjustEdgesLogical & (nearestN != gmrfOrder)) |
			(adjustEdgesLogical & is.null(attributes(N)$cells))
			) {
		N = NNmat(theraster, nearest = gmrfOrder, adjustEdges = adjustEdgesLogical)	
	}
	
	
	paramInfo = list(
			original=param,
			theo = NULL,
			optimal = NULL
	)
	
	
	midcellCoord = c(round(Nx*0.4),round(Ny*0.4)) # the middle cell
	midcell = c(Nx*(Ny-midcellCoord[2]) + midcellCoord[1])
	midVec = sparseMatrix(midcell,1,x=1,
			dims=c(ncol(N),1))
	
	distVecFull =  expand.grid(
			x=seq(1, Nx),
			y=seq(Ny,1))
	buffer = param['shape']+2
	
	distVecFull = t(distVecFull)-midcellCoord 
	distVecFull = sqrt(apply(distVecFull^2,2,sum))	
	
	# if 1 - AR parameter supplied
	if(all(c("oneminusar","shape") %in% names(param))){
		a = 4/(1-param['oneminusar'])	
		
		paramInfo$theo = c(param[c('shape','cellSize','oneminusar')],
				rangeInCells=as.numeric(sqrt(8*param['shape'])/sqrt(a-4)),
				a=as.vector(a), 
				scale = as.numeric(sqrt(a-4))
		)
		if(param['shape'] == 0) {
			paramInfo$theo['rangeInCells'] = as.vector(
					1/paramInfo$theo['scale'])
		}
		
		if ('conditionalVariance' %in% names(param)){
			if(param['shape'] != 0) {
				paramInfo$theo['variance'] =   
						param['conditionalVariance']/
						(pi*param['shape'] *(1-param['oneminusar'])*
							param['oneminusar']^(param['shape'] ))
			} else { # shape zero
				paramInfo$theo['variance'] =  
						4*param['conditionalVariance']/
						((1-param['oneminusar'])*pi)
			}
			paramInfo$theo['conditionalVariance'] = 
					param['conditionalVariance']	
			param['variance'] = paramInfo$theo['variance']
		} else { # conditional variance not supplied
			paramInfo$theo['variance'] = param['variance']
			if(param['shape'] != 0) {
				paramInfo$theo['conditionalVariance']  =   
						paramInfo$theo['variance']*
						(pi*param['shape'] *(1-param['oneminusar'])*
							param['oneminusar']^(param['shape'] ))
			} else {
				paramInfo$theo['conditionalVariance'] =   
						paramInfo$theo['variance']*
						((1-param['oneminusar'])*pi)/4
			}
			
			
		}
		
		
		if(min(c(Nx,Ny)<3*paramInfo$theo['rangeInCells'])){
			warning("grid is ", Nx, " by ", Ny,
					"which may be too small for range ",
					paramInfo$theo['rangeInCells'], " cells, 1- Ar param",
					param['oneminusar'])
		}
		
		paramInfo$theo['range'] = as.numeric(
				paramInfo$theo['rangeInCells']*param['cellSize']
		)
		
		
		#####################	
		# else range supplied
		##########################
	} else if(all(c('range','shape') %in% names(param))){
		
		param['rangeInCells'] = as.numeric(param['range']/param['cellSize'])
		paramInfo$theo = param
		
		
		if(param['shape'] != 0) {
			scale = as.numeric(
					sqrt(8*param['shape'])/param["rangeInCells"])
		} else {
			scale = as.numeric(1/param["rangeInCells"])
		}
		
		a = (scale^2 + 4) 
		paramInfo$theo['oneminusar'] = as.numeric(1-4/a)
		paramInfo$theo['a'] = as.numeric(a)
		
		if(param['shape'] != 0) {
			if ('conditionalVariance' %in% names(param)){
				paramInfo$theo['variance'] =   
						param['conditionalVariance']/
						(pi*param['shape'] *(1-paramInfo$theo['oneminusar'])*
							paramInfo$theo['oneminusar']^(param['shape'] ))
			} else {
				paramInfo$theo['conditionalVariance']  =   
						paramInfo$theo['variance']*
						(pi*param['shape'] *(1-paramInfo$theo['oneminusar'])*
							paramInfo$theo['oneminusar']^(param['shape'] ))
			}
		} else { # shape = 0
			if ('conditionalVariance' %in% names(param)){
				paramInfo$theo['variance'] =  
						4*param['conditionalVariance']/
						((1-paramInfo$theo['oneminusar'])*pi)
			} else {
				paramInfo$theo['conditionalVariance'] =   
						paramInfo$theo['variance']*
						((1-paramInfo$theo['oneminusar'])*pi)/4 
			}
		} # end shape 0
		
		if(min(c(Nx,Ny)<3*param['rangeInCells'])){
			warning("grid is ", Nx, " by ", Ny,
					"which may be too small for range ",
					paramInfo$theo['rangeInCells'], " cells")
		}
		
		
		
		##############
		## not range or oneminusar
		#######################
	} else {
		warning("param must have elements named shape  and either oneminusar or range")
		print(param)
		a = NULL
	}
	
	# build the matrix
	precEntries=getPrec(param['shape'], paramInfo$theo['oneminusar'])
	
	theNNmat = N
	
	theN = (precEntries/paramInfo$theo['conditionalVariance'])[theNNmat@x]
	theN[is.na(theN)] = 0
	theNNmat@x = theN
	
	if(FALSE) {
		bob = solve(forceSymmetric(theNNmat))
		image(matrix(bob[,midcell], ncol(theraster),nrow(theraster)))
	}
	######### optimization to see if there are better matern
	# parameters than the theoretical values
	
	theX = distVecFull * paramInfo$theo['cellSize']
	toKeep = which(theX< 1.75*paramInfo$theo['range'])
	
	if(adjustEdges == FALSE) {	
		varMid = Matrix::solve(Matrix::forceSymmetric(theNNmat),midVec)
		
		ev = data.frame(
				x=theX[toKeep], 
				y=as.vector(varMid[toKeep]))
		
		ev = tapply(ev$y, ev$x, mean)
		paramInfo$empirical = data.frame(x=
						as.numeric(names(ev)),
				y = ev)
		
		paramInfo$empirical = paramInfo$empirical[
				order(paramInfo$empirical$x),
		]
	} else { # adjusting edges
		paramInfo$empirical = data.frame(x = theX[toKeep])	
	}
	# optimal parameters so product of gmrf precision with matern is identity
# don't do it if shape parameter is zero 
	if(paramInfo$theo['shape'] > 0 & paramInfo$theo['shape'] < 4) {
		
		# need stuff on top of file
		nbprec = precEntries[nbn]
		nbprec[is.na(nbprec)] = 0
		otherdist = otherdist * param['cellSize']
		ndist = ndist * param['cellSize']
		
		# inverse only vary shape
		startparam = paramInfo$theo['shape']
		scaleparam = pmax(startparam,0.0001)
		
		objfun = function(qq) {
			qq = c(shape = as.numeric(qq), qq,paramInfo$theo[c('variance','range')])
			
			othermatern = matern(otherdist, param= qq)
			
			nmatern = matern(ndist, param=qq)
			names(nmatern) = names(ndist)
			
			nbprod = as.vector(othermatern %*% nbprec)
			thediag = sum(nmatern[names(precEntries)] * precEntries * numbn[names(precEntries)])
			sum(nbprod^2) + 2*length(othercells)*(thediag-1)^2
		}
		
		optInv = optim(startparam, objfun,
				lower=startparam/4,
				upper=4*startparam,
				control=list(parscale=scaleparam),
				method='L-BFGS-B')
		
		paramInfo$optimalShape = 
				c(optInv$par, 
						paramInfo$theo[c('variance','range',
										'cellSize')])
		
		
		# now estimate variance separately
		objfun = function(qq) {
			qq = c(qq,paramInfo$theo[c('variance','range')])
			othermatern = matern(otherdist, param= qq)
			nbprod = as.vector(othermatern %*% nbprec)
			sum(nbprod^2) 
		}
		
		optInv = optim(startparam, objfun,
				lower=startparam/4,
				upper=4*startparam,
				control=list(parscale=scaleparam),
				method='L-BFGS-B')
		
		nmatern = matern(ndist,
				param=c(optInv$par, 
						paramInfo$theo[c('range','variance')])
		)
		names(nmatern) = names(ndist)
		thediag = sum(nmatern[names(precEntries)] * precEntries * numbn[names(precEntries)])
		
		paramInfo$optimal = c(optInv$par,
				paramInfo$theo[c('cellSize','range')],
				paramInfo$theo['variance']/as.numeric(thediag))
		
		
		paramInfo$empirical$optimal =
				matern(paramInfo$empirical$x,
						param=paramInfo$optimal)
		
		for(D in c('theo', 
				grep("^optimal", names(paramInfo),value=TRUE))
				) {
			
			paramInfo$empirical[,D] =
					matern(paramInfo$empirical$x,
							param=paramInfo[[D]])
			
		}
		
	} # end if shape>0
	
	# edge correction
	
	
	if(adjustEdges!=FALSE){
		if(is.logical(adjustEdges)) {
			adjustEdges='theo'
		}
		
		
		innerCells = attributes(N)$cells$inner
		outerCells = attributes(N)$cells$outer
		
		InnerPrecision = theNNmat[innerCells, innerCells]
		cholInnerPrec =Cholesky(InnerPrecision,LDL=FALSE,perm=TRUE)
		
		#A = x[allCells,-allCells]
		#InnerPrecInvChol = Matrix::solve(Matrix::chol(InnerPrecision))
		#Aic = A %*% InnerPrecInvChol
		# AQinvA = Aic %*% t(Aic)
		
		A = theNNmat[innerCells,outerCells]
		
		Aic = solve(cholInnerPrec, 
				solve(cholInnerPrec,A,system='P'),
				system='L')
		
		paramForM = paramInfo[[adjustEdges]]
		outerCoordsCartesian = xyFromCell(
				theraster, outerCells, spatial=TRUE
		)
		
		precOuter = .Call(
			C_gmrfEdge, 
				as.matrix(Aic),
				outerCoordsCartesian, 
				fillParam(paramForM)[c(
								'range','shape','variance',
								'anisoRatio','anisoAngleRadians','nugget')
				])
		
		theNNmat[outerCells,outerCells] = precOuter
		#	theNNmat = forceSymmetric(theNNmat)
		
	} # end edge adjustment 
	
	theNNmat = Matrix::forceSymmetric(theNNmat)
	
	paramInfo$adjustEdges=adjustEdges
	
	paramInfo$precisionEntries = precEntries
	
#	theNNmat = drop0(theNNmat)
	
	
	attributes(theNNmat)$param =	paramInfo$theo[
			!names(paramInfo$theo) %in% c('a','scale','rangeInCells', 'nugget')
	]	
	attributes(theNNmat)$raster = theraster
	attributes(theNNmat)$info = paramInfo
	
	return(theNNmat)
}

