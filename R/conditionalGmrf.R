conditionalGmrf = function(param,
                           Yvec,Xmat, NN,
                           template=NULL, 
                           mc.cores=1, ...) {
  
  names(param) = gsub("sigmasq","variance",names(param))
  names(param) = gsub("tausq","nugget",names(param))
  

  fixed=as.vector(Xmat %*% 
                    param[paste(colnames(Xmat), "betaHat", sep=".")]
  )
  
  
  Q = maternGmrfPrec(NN,
                     param=c(param[c('oneminusar','shape')], variance=1),
                     ...)
  Qchol = Cholesky(Q, LDL=FALSE)
  # Q/sigmasq = var(U)^(-1)
  
  cholQptau = update(Qchol,parent=Q,mult=1/param['nugget'])
  # Qptau = tausq (1/tausq + Q)
  
  
  #	LofQ = expand(Qchol)$L
  #	pRev = as(ncol(LofQ):1, "pMatrix")		
  #	lQLL =  as( t(LofQ %*% pRev) %*% pRev,'dtCMatrix')
  #	QLL = forceSymmetric(tcrossprod(lQLL))
  
  # QLL = P' L  L' P 
  # QLLorig = Prev P' L   L' P Prev'
  # tausq L^(-1) IcQ L^(-1 ') = var(Y)
  
  #	cholIcQ = Cholesky(QLL,LDL=FALSE,perm=TRUE,
  #			Imult=param['variance']/param['nugget'])
  # need to multiply by param['nugget']
  #	ptwice =   as(expand(cholQLL)$P,'sparseMatrix') %*% 
  #			t(as(pRev,'sparseMatrix')) %*% 
  #			as(expand(Qchol)$P,'sparseMatrix')
  #	ptwice2 = as(ptwice, 'pMatrix')
  #	cholIcQ@perm = as.integer(ptwice2@perm-1)
  
  
  residsOrig = Yvec -	fixed
  EUY = 	as.vector(
    solve(cholQptau, residsOrig,system='A')
  )/param['nugget']
  
  
  
  theidm=c(length(EUY),1)
  varOneCell = function(D) {
    thisD = sparseMatrix(D,1,x=1,dims=theidm)
    solveQ = as.vector(solve(Qchol, thisD))
    #		solveQp = solve(cholIcQ, solveQ,system='P')
    #		thisDp = solve(cholIcQ, thisD ,system='P')
    param['variance_optimal']*(
      solveQ[D] - 
        sum(
          as.vector(solve(cholQptau, thisD)) * 
            solveQ
        )/param['nugget']
    )
  }
  
  
  if(mc.cores==1) {
    thediag = mapply(varOneCell, 1:length(residsOrig))
  } else {
    thediag = parallel::mcmapply(varOneCell, 1:length(EUY),
                                 mc.cores=mc.cores,SIMPLIFY=TRUE)
  }
  
  
  
  VUY = as.vector(
    thediag
  )
  
  
  result=cbind(random=EUY,krigeSd= sqrt(VUY),
               fixed=fixed,predict=fixed+EUY,resids=residsOrig)
  if(!is.null(template)){
    resRast = raster::brick(raster(template), nl=ncol(result))
    names(resRast) = colnames(result)
    values(resRast) = as.vector(result)
    result = resRast		
    
  }
  result
}

