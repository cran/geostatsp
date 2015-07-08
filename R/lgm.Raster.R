setMethod("lgm", 
    signature("formula", "data.frame", "Raster", "data.frame"), 
    function(formula,
        data, grid,
				covariates=NULL,
        shape=1,boxcox=1,nugget=0,
				  expPred=FALSE, nuggetInPrediction=TRUE,
          reml = TRUE, mc.cores=1,
          fixBoxcox=TRUE,
          fixNugget = FALSE,
          ...)
      {
  NN=NNmat(grid)
  csiz = xres(grid)
        
        
  Xmat = model.matrix(formula, data)		  

	if(nrow(Xmat) != ncell(grid))
		warning("dimensions of data and grid are not compatible")
	
  Yvar = all.vars(formula)[1]
  allYvar = grep(paste("^", Yvar, "[[:digit:]]*$",sep=""), names(data), value=TRUE)
  
  Yvec = as.matrix(data[,allYvar, drop=FALSE], drop=FALSE)

  oneminusar = NULL
  if(!fixNugget)
    nugget = NULL
  
  thel = loglikGmrf(Yvec=Yvec,Xmat=Xmat,
                    NN=NN,
                    propNugget=nugget,
                    boxcox=boxcox, fixBoxcox=fixBoxcox,
                    shape=shape,mc.cores=mc.cores,
                    reml=reml,...)
  mle = thel$mle             
  lArray = thel$mlArray
  lMat = thel$mlMat
  

  if (reml){
    chooseLike = 'logLreml'
    m2Like = 'm2logLreml'
  }else{
    chooseLike = 'logLml'
    m2Like = 'm2logLml'
  }
  

  res = list(param = mle)
  if(!is.null(lArray)){
    res$array = lArray
    res$profL = list()
    # nugget
    if(dim(lArray)[2]>1){ # have nugget
      best = apply(lArray[,,m2Like,,drop=FALSE],
          2, which.min) 
      best=arrayInd(best, dim(lArray)[-(2:3)])
      res$profL$propNugget = NULL
      for(D in 1:nrow(best)){
        res$profL$propNugget = rbind(
            res$profL$propNugget,
            lArray[best[1], D, 
                c('propNugget',chooseLike),
                best[2]])
      }
    } # end have nugget
     
    if(dim(lArray)[4]>1){ # have oneminusar
      best = apply(lArray[,,m2Like,,drop=FALSE],
          4, which.min) 
      best=arrayInd(best, dim(lArray)[-(3:4)])
      
      res$profL$oneminusar = NULL
      for(D in 1:nrow(best)){
        res$profL$oneminusar = rbind(
            res$profL$oneminusar,
            lArray[best[1], best[2], 
                c('oneminusar',chooseLike,'range'),
                D])
      }
    } # end have oneminusar
    res$profL$range = res$profL$oneminusar[,c(3,2)]
    res$profL$oneminusar = res$profL$oneminusar[,c(1,2)]
    
    if(all(dim(lArray)[c(2,4)]>1)){ # have both
      x=lArray[1,,'propNugget',1]
      orderx = order(x)
      y = lArray[1,1,'range',]
      ordery=order(y)
      
      res$profL$twoDim = list(
          x=lArray[1,orderx,'propNugget',1],
          y=lArray[1,1,'range',ordery],
          z=apply(
          lArray[,,chooseLike,,drop=FALSE],
          c(2,4), max, na.rm=TRUE)[orderx,ordery],
      oneminusar=lArray[1,1,'oneminusar',ordery]      
    )
 
    } # end have both
  }

  
  res$data = data
  res$model$reml = reml
  res$model$trend = formula
 
  # summary table
  scovariates = gsub(
      'stdErr\\.','',
      grep("^stdErr\\.", names(mle), value=TRUE)
  )
  
  srownames = c('sdNugget','sdSpatial','range','shape')
  scolnames = c("estimate", "stdErr", "ci0.005", "ci0.995", "ci0.025", "ci0.975", 
      "ci0.05", "ci0.95", "ci0.1", "ci0.9", "pval", "Estimated")
  ress = as.data.frame(
      matrix(
          NA,
          length(scovariates) + length(srownames),
          length(scolnames),
          dimnames = list(
              c(scovariates, srownames),
              scolnames
              )
          )
      )
   ress[c('sdNugget','sdSpatial','range','shape'),'estimate'] = 
       c(sqrt(mle[c('nugget','variance')]),
        mle[c('range','optimalShape.shape')]   
       )
   ress[c('sdNugget','sdSpatial','range','shape'),'Estimated'] =
       c(
           fixNugget, TRUE, TRUE, FALSE
           
           )
   ress[scovariates,'Estimated']  = TRUE   
   ress[scovariates,'estimate']  = mle[scovariates]   
   ress[scovariates,'stdErr']  = mle[paste("stdErr.",scovariates,sep="")]   
   
       
  
  res$summary = ress

	
  return(res)
  

}
)

