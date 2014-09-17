

loglikGmrfGivenQ = function(
		propNugget,
  Ry, Y, 
  Rx, X,
  Q, Qchol=NULL, 
  detQ=NULL) {
  
  
  # propNugget =    tau^2/xi^2
  
	if(is.vector(Y))	
		Y = as.matrix(Y)
	
  N= nrow(Y) + c(ml=0,reml=-ncol(X))

  if(is.null(Qchol)) {
	  Qchol = Cholesky(Q, LDL=FALSE)
  }
  if(is.null(detQ)) {
	  detQ= 2*determinant(Qchol,logarithm=TRUE)$modulus
  }
  
  if(propNugget>0){

	  xisqtausq = 1/propNugget
  
	  Vchol = update(Qchol, Q, mult=xisqtausq)
	  
	  Vy = solve(Vchol, Y)
	  Vx = solve(Vchol, X)
	  
    
    XprecXinv = solve(Matrix::crossprod(Rx,Vx))
    
    betaHat =  
      XprecXinv %*% 
        Matrix::crossprod(Rx,Vy) 
    
    rownames(betaHat) = colnames(X)

	R = diag(
			Matrix::crossprod(Ry,Vy) - 
					Matrix::crossprod(Ry, Vx) %*% betaHat
	)
 
	
	constHat=tausq = outer(N,R,"/")
	xisq = tausq/propNugget
	
    logDetVar = as.numeric(
      2*determinant(Vchol,logarithm=TRUE)$modulus
	)
        
  } else { # no nugget 
    
	  
    XprecX = Matrix::crossprod(X,Rx)
    XprecXinv = solve(XprecX)
    
    betaHat = XprecXinv %*% 
                           Matrix::crossprod(X , Ry)
    rownames(betaHat) = colnames(X)
    
    resids  =  Y - X %*% betaHat 
    R =  t(Matrix::crossprod(resids,Q)) * resids
	R = apply(R, 2,sum)
	
    logDetVar=0
	
    constHat = xisq = t(outer(R,N,"/"))
	tausq=0*constHat
	
    # matrix operations
    # Qchol, cholCovMat = Matrix::chol(covMat)
    # XL, cholCovInvX = Matrix::solve(cholCovMat, covariates)
    # XprecX, cholCovInvXcross = Matrix::crossprod(cholCovInvX)
    # XprecXinv, cholCovInvXcrossInv = Matrix::solve(cholCovInvXcross)
  }
  
  logDetVar = logDetVar - detQ + N - N*log(N)
  m2logL = logDetVar + outer(N,log(R), "*")
  
  
  variances = abind::abind(tausq = tausq, xisq = xisq,
		  along=3)
  
  parInfo = attributes(Q)$param
  for(D in c('theo', 
		 grep("^optimal", names(parInfo),value=TRUE))) {

	 sigmasq = xisq * parInfo[[D]]['variance']
	 variances = abind::abind(variances, sigmasq, along=3)
	 dimnames(variances)[[3]][dim(variances)[3]] = 
			 paste('sigmasq_',D,sep='')
 }
  
  variances = apply(variances, 2, function(x) {
			  res = as.vector(x)
			  names(res)=t(outer(
					  colnames(x),rownames(x),
					  paste,sep='.'))
			  res})
  

  
  
  
  m2logL['reml',] =	m2logL['reml',] +  determinant(
    XprecXinv,logarithm=TRUE)$modulus
  
  rownames(m2logL) = paste("m2logL.",
		  rownames(m2logL),sep='')
  
  logL = -m2logL/2
  rownames(logL) = gsub("^m2","",rownames(logL))
  
  sebeta = diag(XprecXinv)
  sebeta = array(sebeta, 
		  c(length(sebeta),dim(constHat)),		
		  dimnames=c(
				  list(rownames(betaHat)),
				  dimnames(constHat)))
  
  constHat2 = array(constHat, dim(sebeta)[c(2,3,1)],
		  dimnames=dimnames(sebeta)[c(2,3,1)])
  constHat2 = aperm(constHat2, c(3,1,2))
  sebeta = sqrt(sebeta *constHat2)
  names(sebeta)= paste('se',rep(names(betaHat),2), 
                       names(sebeta),sep='.')

  sebeta = apply(sebeta, 3, function(x) {res = as.vector(x);names(res)=outer(rownames(x),colnames(x),paste,sep='.se.');res})
  									   
#  varbetahat = abind::abind(
#    as.matrix(XprecXinv)/constHat[1],
#    as.matrix(XprecXinv)/constHat[2],
#    along=3)
  
#  dimnames(varbetahat) = 
#    list(names(betaHat),names(betaHat),
#         names(constHat))
  

  betaHat = as.matrix(betaHat)
  rownames(betaHat) = paste(rownames(betaHat), "betaHat", sep=".")

   result =rbind(m2logL, logL,
            variances,
            as.matrix(betaHat), 
            sebeta,
			propNugget=propNugget)#,
#            logDetVar=as.numeric(logDetVar))
  
  
#  attributes(result)$varbetahat = varbetahat
  
  result
  
}

loglikGmrfOneRange = function(
  oneminusar,
  Yvec, Xmat, NN, propNugget=0,
  shape=1,  
  adjustEdges=FALSE) {
    
    Q =  maternGmrfPrec(NN,
                         param=c(shape=as.vector(shape),
                                 oneminusar=as.vector(oneminusar[1]),
                                 conditionalVariance=1),
		adjustEdges=adjustEdges)
    
	
  
	Qchol = Cholesky(Q,LDL=FALSE)
 
  
  detQ = 2*as.numeric(determinant(Qchol,
                                    logarithm=TRUE)$modulus)
  
  Rx = Q %*% Xmat
  Ry = Q %*% as.matrix(Yvec)
  
  argList = list(
  	Ry=Ry, Y=as.matrix(Yvec), 
  	Rx=Rx, X=Xmat,
 	Q=Q, Qchol=Qchol, 
  detQ=detQ)  
  
thepar=attributes(Q)$param
   
  if(length(propNugget)==1) {
    
    argList$propNugget=propNugget
    
    res = do.call(loglikGmrfGivenQ,argList)
    
    res = as.matrix(res)
    
  } else {
    
    res = mapply(loglikGmrfGivenQ,
                 propNugget=propNugget,			
                 MoreArgs=argList,SIMPLIFY=FALSE
    )
    names(res) = paste("propNugget=", propNugget,sep="")
	res = simplify2array(res)    

}	  




res = abind::abind(res, 
		oneminusar=array(oneminusar, dim=dim(res)[-1]), 
		range = array(thepar$theo['range'],dim=dim(res)[-1]), 
		shape = array(thepar$theo['shape'],dim=dim(res)[-1]), 
		shape_optimal = array(thepar$optimal['shape'],dim=dim(res)[-1]), 
		along=1)

res = drop(res)

attributes(res)$Qinfo = thepar

res

}

#  for(D in c("theo", "optimal")#,"optimalSameShape",
	if(F
#		  "optimalInverse", "optimalInverseSameShape",
#		  "optimalInverseSameRange", 'optimalWithNugget',
#		  'optimalInverseShape', 'optimalShape')
  ) {
	  
  
	newres = rbind(range=
					thepar[[D]]['range'],
			shape=thepar[[D]]['shape'],
			sigmasq.ml = res['xisq.ml',] * 
					thepar[[D]]['variance'],
			sigmasq.reml = res['xisq.reml',] * 
					thepar[[D]]['variance']
	)

	rownames(newres) = paste(rownames(newres),
			D, sep="_")	
	  
	res = rbind(res, newres)
	  
  }
  
if(F){	res = rbind(res,
		  tausq.ml_optimalWithNugget = res['tausq.ml',] + 
				  thepar$optimalWithNugget['nugget'],
		  tausq.reml_optimalWithNugget = 
				  res['tausq.reml',] + 
				  thepar$optimalWithNugget['nugget']
		  )
  
	rownames(res) = gsub("_theo$", "", 
			rownames(res))		  
		  
}




loglikGmrf = function(
  oneminusar,
  Yvec, Xmat, NN, propNugget=0,
  shape=1, 
  adjustEdges=FALSE,
  mc.cores=1) {
  
  
  argList = list(Yvec=as.matrix(Yvec), 
                 Xmat=as.matrix(Xmat), 
				 NN=NN,
                 propNugget=propNugget,
                 shape=shape,
                 adjustEdges=adjustEdges
  )
  
  if(mc.cores>1) {
    myapply= function(...){
      parallel::mcmapply(...,
                         mc.cores=mc.cores)
    }
  } else {
    myapply= mapply
  }
  
    res=myapply(loglikGmrfOneRange,
                oneminusar=oneminusar,
                MoreArgs=argList,
                SIMPLIFY=FALSE
    )
    names(res) = 
      paste("oneminusar=",oneminusar,sep="")
  
  if(class(res[[1]])!= 'character'){
    Qinfo = attributes(res[[1]])$Qinfo
    
    res= simplify2array(res)
    
    attributes(res)$Qinfo = Qinfo
    
  }
  res
}



summaryGmrfFit= function(x) {
  UseMethod("summaryGmrfFit")  
}

summaryGmrfFit.array = function(x) {
  
  
  x2 = aperm(x,c(3,2,1))
  x2 = matrix(c(x2), ncol=dim(x2)[3])
  colnames(x2) = dimnames(x)[[1]]
  res = summaryGmrfFit.matrix(x2,npar=2)
  return(res)
}


summaryGmrfFit.matrix = function(x,npar=1) {
  
  if(any(rownames(x)=='logL.ml')){
    x = t(x)
  }
  result = list()
  someL = c('ml','reml')
  for(D1 in someL) {
    D=paste('logL.',D1,sep='')
    MLErow = x[which.max(x[ ,D]),]  
    withinCI = which(-2*(x[,D] - MLErow[D]) < qchisq(0.95, 2))
    
    xsub = x[,grep("^se\\.",colnames(x),invert=TRUE),drop=FALSE]
    parCI= t(apply(
      xsub[withinCI, ,drop=FALSE],
      2,range))
    colnames(parCI) = c('q0.025', 'q0.975')
    
    parMat = cbind(mle=MLErow[colnames(xsub)],parCI,se=NA) 
    notD = someL[someL != D1]
    notD = grep(
      paste("\\.",notD, "$",sep=''),
      rownames(parMat), 
      invert=TRUE,value=TRUE)
    
    parMat = parMat[notD,]
    rownames(parMat) = gsub(paste("\\.", D1, "$",sep=""),
                            "",rownames(parMat))
    
    se = grep("\\.se$", rownames(parMat), value=TRUE)
    seInt = parMat[grep('^\\(', se, value = T),1]
    seX = parMat[grep('^x', se, value = T),1]
    parMat['(Intercept).betaHat',4] = seInt
    parMat['x.betaHat',4] = seX
    parMat['x.betaHat',2] = as.numeric(parMat['x.betaHat',1]-2*seX)
    parMat['x.betaHat',3] = as.numeric(parMat['x.betaHat',1]+2*seX)
    parMat['(Intercept).betaHat',2] = as.numeric(parMat['(Intercept).betaHat',1]-2*seInt)
    parMat['(Intercept).betaHat',3] = as.numeric(parMat['(Intercept).betaHat',1]+2*seInt)
    parRid = parMat[which(rownames(parMat) %in% se),]
    judgRid = !(rownames(parMat) %in% rownames(parRid))
    parKeep = parMat[which(judgRid == T),]
    dele = c(grep("^shape",rownames(parKeep)), grep("^sigmasq",rownames(parKeep)))
    parKeep[dele,2] = NA
    parKeep[dele,3] = NA
    result[[D1]] = parKeep
    
  }
  return(result)		
}

