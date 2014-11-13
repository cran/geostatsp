# 
# setMethod("lgm", 
# 		signature("formula", "data.frame", "Raster", "data.frame"), 
lgm.Raster=	function(formula,data, grid,
				covariates=NULL,
                  shape=1,boxcox=1,nugget=0,
				  expPred=FALSE, nuggetInPrediction=TRUE,
				  reml = TRUE,
				  mc.cores=1,		
				  oneminusar=seq(0.01, 0.4,len=4), 
                  range=NULL,  
                  ...) {

	Xmat = model.matrix(formula, data)		  

	if(nrow(Xmat) != ncell(grid))
		warning("dimensions of data and grid are not compatible")
	
  NN=NNmat(grid)
  csiz = xres(grid)

  Yvar = all.vars(formula)[1]
  allYvar = grep(paste("^", Yvar, "[[:digit:]]*$",sep=""), names(data), value=TRUE)
  
  Yvec = as.matrix(data[,allYvar, drop=FALSE], drop=FALSE)
  reXmat = Xmat
#	colnames(reXmat) = c('(Intercept)', 'x') 
  thel = loglikGmrf(oneminusar=oneminusar,
                    Yvec=Yvec,Xmat=reXmat,
                    NN=NN,propNugget=nugget,
                    shape=shape,mc.cores=mc.cores,...)
  thesummary = list()
  if (reml){
    chooseLike = 'logL.reml'
    m2Like = 'm2logL.reml'
  }else{
    chooseLike = 'logL.ml'
    m2Like = 'm2logL.ml'
  }
  


  if (all(nugget == 0)) {
    thesummary$profL$propNugget = 0
    propNug = 0
    omar = thel['oneminusar',]
    rangeValue = csiz*sqrt(2*shape*(1-omar)/omar)
    rangeInCellsValue = rangeValue/csiz
    rangeLike = thel[chooseLike,]
    rangem2Like = thel[m2Like,]
    rangePack = cbind(omar, rangeValue, rangeInCellsValue, rangeLike, rangem2Like)
    rownames(rangePack) = NULL
    colnames(rangePack) = c('oneminusar', 'range', 'rangeInCells',chooseLike, m2Like)
    
    thesummary$profL$range = as.data.frame(rangePack)
  }
  else{
    # $propNugget
		print(dimnames(thel))
    propNug = thel['propNugget',,1]
    propLike =  thel[chooseLike,,1]
    propm2Like =  thel[m2Like,,1]
    propPack = cbind(propNug, propLike, propm2Like)
    rownames(propPack) = NULL
    colnames(propPack) = c('propNugget', chooseLike, m2Like)
    
    thesummary$profL$propNugget = as.data.frame(propPack)
    
    
    #$range
    omar = thel['oneminusar',1,]
    rangeValue = csiz*sqrt(2*shape*(1-omar)/omar)
    rangeInCellsValue = rangeValue/csiz
    rangeLike = thel[chooseLike,1,]
    rangem2Like = thel[m2Like,1,]
    rangePack = cbind(omar, rangeValue, rangeInCellsValue, rangeLike, rangem2Like)
    rownames(rangePack) = NULL
    colnames(rangePack) = c('oneminusar', 'range', 'rangeInCells',chooseLike, m2Like)
    
    thesummary$profL$range = as.data.frame(rangePack)
  }
  #$twoDim
  thesummary$profL$twoDim = list()
  thesummary$profL$twoDim$oneminusar = omar
  thesummary$profL$twoDim$range = rangeValue
  thesummary$profL$twoDim$rangeInCells = rangeInCellsValue
  thesummary$profL$twoDim$propNugget = propNug
  thesummary$profL$twoDim$array = thel
  
 
  thesummary$data = data
  thesummary$model$reml = reml
  thesummary$model$trend = formula
 
  thesummary$summary = summaryGmrfFit(thel)[,,c('ml','reml')[reml+1]]

	
  thesummary$param = thesummary$summary[,'mle']
  
  mleparam = thesummary$param
  

  
  if(mleparam['propNugget']>0) {
    thesummary$predict = conditionalGmrf(
      param=mleparam,
      Yvec=Yvec,Xmat=Xmat,
      template=grid, NN=NN,
      mc.cores=mc.cores,...)
  } else {
    thesummary$predict = raster::brick(
      raster(grid), nl=ncol(Yvec)+1)
    names(thesummary$predict) = c('fixed',paste('random', colnames(Yvec), sep="."))
	values(thesummary$predict) = NA
    values(thesummary$predict[['fixed']]) =
      Xmat %*% mleparam[
        paste(colnames(reXmat),".betaHat",sep='')]
	for(D in colnames(Yvec)) {
    values(thesummary$predict[[paste('random', D, sep=".")]]) =
      Yvec[,D] - values(thesummary$predict[['fixed']])
	}
  }
  if(nlayers(thesummary$predict)==2)
	  names(thesummary$predict) = c("fixed","random")
  
  return(thesummary)
}


setMethod("lgm", 
		signature("formula", "data.frame", "Raster", "data.frame"), 
		lgm.Raster)



