loglikGmrfGivenQobjFun = function(
    propNugget,
    reml,
    ...
    ) {

      resFull = loglikGmrfGivenQ(
          propNugget=propNugget,
          reml=reml, 
          ...
              )
       res = resFull[
           paste('logL.',c('ml','reml')[1+reml], sep=''),
           1]
       attributes(res)$full = resFull
       res
    }

loglikGmrfObjFunBc = function(
    onebc, 
    Y, X,
    propNugget,
    Qchol, Vchol, 
    Vx,
    cholXXLxLx,
    boxcoxInterval=NULL, 
    sumLogY=NULL,
    reml=TRUE
) {
      if(abs(onebc)<0.001) {
        Ybc = log(Y) 
      } else  { #boxcox far from 0 and 1
        Ybc <- ((Y^onebc) - 1)/ onebc 
      }
      if(is.null(sumLogY))
        sumLogY = sum(log(Y[,1]))

      twoLogJacobian = as.numeric(2*(onebc-1)*sumLogY) 
      
      if(propNugget>0){
        YprimeY = apply(Y^2,2,sum)
        YprimeX = crossprod(Y,X)
        Vy = solve(Vchol, Y,system='L')
        XYLxLy = t(YprimeX) - crossprod(Vx,Vy)/propNugget 
      } else { # no nugget 
        Vy = expand(Qchol)$L %*% Y
        XYLxLy = crossprod(Vx,Vy)
        YprimeY = XprimeX = 0
      }
      
      Lc = solve(cholXXLxLx, XYLxLy)
      LcLc = apply(Lc^2, 2, sum)
      VyVy = apply(Vy^2,2,sum)
      
      if(propNugget>0){
        ssq = YprimeY - VyVy/propNugget  - LcLc
      } else { # no nugget 
        ssq = VyVy  - LcLc
      }
      
      N = nrow(Y)
      if(reml) {
        N = N - ncol(X)
      }
      
      # return logL + constant, not minus2 log L
      N * log(ssq[1]) - twoLogJacobian
}
    
loglikGmrfGivenQ = function(
        propNugget,
        Y,X,
        Q, 
        Qchol=NULL, 
        detQ=NULL,
        boxcoxInterval=NULL,
        reml=TRUE
    ) {
      # propNugget =    tau^2/xi^2
      # if Qchol is specified, 
      # it is assumed that Y and X have already
      # been rotated
      # reml is only used for choosing box-cox parameter
      if(is.vector(Y))	
        Y = as.matrix(Y)
      
      if(is.null(Qchol)) {
        Qchol = Cholesky(Q, LDL=FALSE)
        Y=solve(Qchol, Y, system='P')
        X=solve(Qchol, X, system='P')
      }
      if(is.null(detQ)) {
        detQ= 2*determinant(Qchol,logarithm=TRUE)$modulus
      }
      
      if(propNugget>0){
        XprimeX = crossprod(X)
        
        xisqtausq = 1/propNugget
        Vchol = update(Qchol, Q, mult=xisqtausq)
        detLhalf = determinant(Vchol,logarithm=TRUE)$modulus
        Vx = solve(Vchol, X,system='L')
        
        XXLxLx = XprimeX - xisqtausq * crossprod(Vx) 
      } else { # no nugget 
        Vchol = NULL
        detLhalf = 0
        Vx = expand(Qchol)$L %*% X
        XXLxLx = crossprod(Vx)
      }
      cholXXLxLx = t(chol(XXLxLx))
      detXXhalf = Matrix::determinant(
          cholXXLxLx,logarithm=TRUE)$modulus
      
      if(length(boxcoxInterval)==2){ # run an optimizer
        
        boxcox = rep(NA, ncol(Y))
        sumLogY = apply(log(Y), 2, sum)
        
        for(Dvar in 1:ncol(Y)) {
          
          boxcox[Dvar] =  optimize(
              loglikGmrfObjFunBc, 
              range(boxcox),
              maximum=FALSE, 
              Y=Y[,Dvar], 
              X=X,
              propNugget=propNugget,
              Qchol=Qchol,
              Vchol=Vchol, 
              Vx=Vx,
              cholXXLxLx=cholXXLxLx,
              reml=reml,
              sumLogY=sumLogY[Dvar]
          )$minimum
          #found boxcox
          if(abs(boxcox[Dvar])<0.001) {
            Y[,Dvar] = log(Y[,Dvar]) 
          } else  { #boxcox far from 0 and 1
            Y[,Dvar] <- ((Y[,Dvar]^boxcox[Dvar]) - 1)/ boxcox[Dvar] 
          }
        }
        # Y is now box-cox transfomed
        twoLogJacobian = as.numeric(2*(boxcox-1)*sumLogY) 
      } else { # no boxcox
        boxcox=rep(1, ncol(Y))
        twoLogJacobian=rep(0, ncol(Y))
      }
      
      if(propNugget>0){
        YprimeY = apply(Y^2,2,sum)
        YprimeX = crossprod(Y,X)
        
        Vy = solve(Vchol, Y,system='L')
        XYLxLy = t(YprimeX) - crossprod(Vx,Vy)/propNugget 
      } else { # no nugget 
        Vy = expand(Qchol)$L %*% Y
        XYLxLy = crossprod(Vx,Vy)
      }

      Lc = solve(cholXXLxLx, XYLxLy)
      LcLc = apply(Lc^2, 2, sum)
      VyVy = apply(Vy^2,2,sum)
      
      if(propNugget>0){
        ssq = YprimeY - VyVy/propNugget  - LcLc
      } else { # no nugget 
        ssq = VyVy  - LcLc
      }
      
      # prepare the result
      N = nrow(Y)
      Ncov = ncol(X)
      Nreml = N - Ncov
      result = list(
          m2logL.ml = N * log(ssq)+ 2*detLhalf - detQ + 
              N - N*log(N) - twoLogJacobian, 
          m2logL.reml = Nreml * log(ssq)+ 2*detLhalf - 
              detQ -2*detXXhalf + Nreml - 
              Nreml * log(Nreml) - twoLogJacobian
      )
      
      result$logL.ml = -result$m2logL.ml/2
      result$logL.reml = -result$m2logL.reml/2
      
      const.ml = ssq/N
      const.reml = ssq/Nreml
      
      
      if(propNugget>0){
        result$tausq.ml = const.ml
        result$tausq.reml = const.reml
        result$xisq.ml = result$tausq.ml/propNugget
        result$xisq.reml = result$tausq.reml/propNugget
      } else { # no nugget 
        result$tausq.ml = result$tausq.reml = 0
        result$xisq.ml = const.ml
        result$xisq.reml =const.reml
      }
      parInfo = attributes(Q)$param
      result$sigmasq_theo.ml = result$xisq.ml * parInfo$theo[['variance']]  
      result$sigmasq_theo.reml = result$xisq.reml * parInfo$theo[['variance']]  
      
      result$sigmasq_optimalShape.ml = 
          result$xisq.ml * parInfo$optimalShape[['variance']]  
      result$sigmasq_optimalShape.reml = 
          result$xisq.reml * parInfo$optimalShape[['variance']]  
      
      result$betaHat = as.matrix(solve(t(cholXXLxLx),Lc))
      rownames(result$betaHat ) = paste(
          rownames(result$betaHat), ".betaHat",sep='')
      
      seBetaHat = diag(chol2inv(cholXXLxLx))
      seBetaHat = matrix(seBetaHat, nrow=length(seBetaHat),
          ncol=ncol(Y))
      rownames(seBetaHat) = paste(colnames(X), '.se',sep='')
      
      result$seBetaHat.ml = seBetaHat * matrix(
          const.ml,
          nrow(seBetaHat), ncol(seBetaHat), 
          byrow=TRUE)
      rownames(result$seBetaHat.ml) = paste(
          rownames(result$seBetaHat.ml), '.ml',sep=''
      )
      
      
      result$seBetaHat.reml = seBetaHat * matrix(
          const.reml,
          nrow(seBetaHat), ncol(seBetaHat), 
          byrow=TRUE)
      rownames(result$seBetaHat.reml) = paste(
          rownames(result$seBetaHat.reml), '.reml',sep=''
      )
      
      result$propNugget = rep(propNugget, ncol(Y))
      result$boxcox = boxcox
      
      result  = do.call(rbind, result)
      
      result
}    
    
loglikGmrfGivenQold = function(
	propNugget,
  Ry, Y, 
  Rx, X,
  Q, Qchol=NULL, 
  detQ=NULL, 
  boxcoxInterval=NULL, 
  reml=TRUE,
  sumLogY = NULL) {
  
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
	  
	Vx = solve(Vchol, X)
   	XprecXinv = solve(Matrix::crossprod(Rx,Vx))

	if(length(boxcoxInterval)==2){ # run an optimizer

		boxcox = rep(NA, ncol(Y))
		for(Dvar in 1:ncol(Y)) {
			bchere = optimize(oneL, range(boxcox),
					maximum=TRUE, Y=Y[,Dvar], sumLogY=sumLogY[Dvar])$maximum
			boxcox[Dvar] = bchere
		 #found boxcox
			if(abs(bchere)<0.001) {
				Y[,Dvar] = log(Y[,Dvar]) 
			} else  { #boxcox far from 0 and 1
				Y[,Dvar] <- ((Y[,Dvar]^bchere) - 1)/ bchere 
			}
		}
		Ry = Q %*% Y
		# Y is now box-cox transfomed
		twoLogJacobian = as.numeric(2*(boxcox-1)*sumLogY) 
	} else { # no boxcox
		boxcox=NULL
		twoLogJacobian=0
	}
	Vy = solve(Vchol, Y)
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
	  
  	if(length(boxcoxInterval)==2){ # run an optimizer
	  oneL = function(onebc, Y, sumLogY) {
		  if(abs(onebc)<0.001) {
				  Ybc = log(Y) 
		  } else  { #boxcox far from 0 and 1
				  Ybc <- ((Y^onebc) - 1)/ onebc 
		  }
			  twoLogJacobian = as.numeric(2*(onebc-1)*sumLogY) 
			  
			  Rybc = Q %*% Ybc
			  
			  betaHat = XprecXinv %*% 
					  Matrix::crossprod(X , Rybc)
			  
			  resids  =  Ybc - X %*% betaHat 
			  R =  t(Matrix::crossprod(resids,Q)) * resids
			  R = apply(R, 2,sum)
			  
			  logDetVar = - detQ + N - N*log(N)
			  m2logL = logDetVar + outer(N,log(R), "*") - twoLogJacobian
			  m2logL[ c('ml','reml')[1+reml] ]
	  }
	boxcox = rep(NA, ncol(Y))
	for(Dvar in 1:ncol(Y)) {
		bchere = optimize(oneL, range(boxcox),
				maximum=TRUE, Y=Y[,Dvar], sumLogY=sumLogY[Dvar])$maximum
		boxcox[Dvar] = bchere
		  #found boxcox
		  if(abs(bchere)<0.001) {
			  Y[,Dvar] = log(Y[,Dvar]) 
		  } else  { #boxcox far from 0 and 1
			  Y[,Dvar] <- ((Y[,Dvar]^bchere) - 1)/ bchere 
		  }
	}
	Ry = Q %*% Y
	twoLogJacobian = as.numeric(2*(boxcox-1)*sumLogY) 	
  } else {
	  boxcox=NULL
	  
	twoLogJacobian=0  
  } # end box cox
  
	  
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
  } # end no nugget
  
  logDetVar = logDetVar - detQ + N - N*log(N)
  m2logL = outer(N,log(R), "*")
  m2logL =  m2logL  + (logDetVar  - c(twoLogJacobian,twoLogJacobian))
  
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

  sebeta = apply(sebeta, 3, function(x) {
        res = as.vector(x)
        names(res)=outer(rownames(x),colnames(x),paste,sep='.se.')
        res}
  )
  									   
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
			propNugget=propNugget, boxcox=rep(boxcox,2))#,
#            logDetVar=as.numeric(logDetVar))
  
  
#  attributes(result)$varbetahat = varbetahat
  
  result
  
}

loglikGmrfOneRange = function(
  oneminusar,
  Yvec, Xmat, NN, propNugget=0,
  shape=1,  boxcoxInterval=NULL,
  reml=TRUE,
  sumLogY = NULL,
  adjustEdges=FALSE,
  optimizer=FALSE) {
    
    Q =  maternGmrfPrec(NN,
                         param=c(shape=as.vector(shape),
                                 oneminusar=as.vector(oneminusar[1]),
                                 conditionalVariance=1),
		adjustEdges=adjustEdges)

  thepar=attributes(Q)$param

if(FALSE){
  # eventually, this will be done in C
  xisqtausq = 1/propNugget
  Nnugget=length(xisqtausq)
  Nobs = nrow(Xmat)
  Ncov = ncol(Xat)
  if(is.vector(Yvec))	
    Yvec = as.matrix(Yvec)
  Ny = ncol(Yvec)
  resC=.C(
      'logLikGmrfGivenQ',
      xisqtausq, Nnugget,
      Yvec, as.matrix(Xmat), 
      Nobs, Ny, Ncov,
      Q,
      betaHat = rep(0.0, Ncov*Ny*Nnugget),
      cholXVX = rep(0.0, Ncov*Ncov*Nnugget),
      determinants = rep(0.0,3*Nnugget),
      ssq = rep(0.0, Ny*Nnugget)
  )
}
  
	Qchol = Cholesky(Q,LDL=FALSE)
  detQ = 2*as.numeric(determinant(Qchol,
                                    logarithm=TRUE)$modulus)
  argList = list(
    Y=solve(Qchol, as.matrix(Yvec), system='P'),
    X=solve(Qchol, as.matrix(Xmat), system='P'),
    Q=Q, Qchol=Qchol, 
    detQ=detQ,
    boxcoxInterval=boxcoxInterval)

  if(optimizer) { # use an optimizer

    
    res = optimize(
        loglikGmrfGivenQobjFun,
        lower=0,upper=100,
        maximum=TRUE,
    # arguments passed to loglikGmrfGivenQ
        Y=argList$Y, 
        X=argList$X,
        Q=argList$Q, 
        Qchol=argList$Qchol, 
        detQ=argList$detQ,
       boxcoxInterval = boxcoxInterval,
        reml=reml#, sumLogY=sumLogY
        )
    res = attributes(res$objective)$full
  } else { # evaluate for a sequence of propNugget
  
  
  if(length(propNugget)==1) {
    
    argList$propNugget=propNugget
    
    res = do.call(loglikGmrfGivenQ,argList)
    
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

}
attributes(res)$Qinfo = thepar

res

}


loglikGmrf = function(
  oneminusar,
  Yvec, Xmat, NN, propNugget=0,
  shape=1, boxcox=1,
  reml=TRUE,
  adjustEdges=FALSE,
  optimizer=FALSE,
  mc.cores=1) {
  	
	if(length(boxcox)>1) {
		boxcoxInterval = range(boxcox)
		sumLogY = apply(as.matrix(Yvec), 2, function(qq) sum(log(qq)))	
		if(is.nan(sumLogY))
			warning("boxcox shouldnt be used with negative data")
	} else { # a single boxcox parameter (perhaps 1)
		boxcoxInterval=NULL
		sumLogY=NULL
		if(abs(boxcox-  1 ) > 0.001) {
			if(abs(boxcox<0.001)) {
				Yvec = log(Yvec) 
			} else  { #boxcox far from 0 and 1
				Yvec <- ((Yvec^boxcox) - 1)/boxcox
			}
		}
	}
	
  	argList = list(Yvec=as.matrix(Yvec), 
                 Xmat=as.matrix(Xmat), 
				 NN=NN,
                 propNugget=propNugget,
                 shape=shape,
                 adjustEdges=adjustEdges,
		 boxcoxInterval=boxcoxInterval, sumLogY=sumLogY,
		 reml=reml,
     optimizer=optimizer
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
  
  
  x2 = aperm(x,seq(length(dim(x)),1))
  x3 = array(x2, c(prod(dim(x2)[1:2]),dim(x2)[-(1:2)]),
		  dimnames=c(
				  list(as.vector(outer(dimnames(x2)[[1]], dimnames(x2)[[2]], paste))),
				  dimnames(x2)[-(1:2)]))
  if(length(dim(x3))==2){
	  res = summaryGmrfFit.matrix(x3,npar=2)
  } else if(length(dim(x3))==3){
	res=NULL
	for(D  in seq(dim(x3)[2],1)){
		res = abind::abind(summaryGmrfFit.matrix(x3[,D,]), 
				res, along=4)
	}
	dimnames(res)[[4]]=dimnames(x3)[[2]]
  }
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
	betahats = gsub("se$", "betaHat", se)
	
	parMat[betahats,"se"] = parMat[se,"mle"]
	

	parMat[betahats,"q0.025"] = as.numeric(parMat[betahats,"mle"]-2*parMat[betahats,"se"])
	parMat[betahats,"q0.975"] = as.numeric(parMat[betahats,"mle"]+2*parMat[betahats,"se"])
	
 
    parKeep = parMat[!rownames(parMat) %in% se,]
	
    dele = grep("^shape|^sigmasq|^xisq|^tausq",rownames(parKeep))
    parKeep[dele,grep("^q0\\.", colnames(parKeep))] = NA
    result[[D1]] = parKeep
    
  }
  result = abind::abind(result[[1]], result[[2]], along=3)
  dimnames(result)[[3]] = someL
  return(result)
}


