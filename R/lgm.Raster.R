
setMethod("lgm", 
		signature("ANY", "Raster"), 
function(formula,data,
				covariates=NULL,
                  shape=1,boxcox=1,nugget=0,
				  newdata=data,
				  expPred=FALSE, nuggetInPrediction=TRUE,
				  reml = TRUE,
				  mc.cores=1,		
				  oneminusar=seq(0.01, 0.4,len=4), 
                  range=NULL, 
                  ...) {
  

			  
			  
  if(!is.null(covariates)){
    covariates = stackRasterList(covariates,template=data)
    data = stack(data,covariates)
  }
  NN=NNmat(data)
  csiz = xres(data)
  data = as.data.frame(data)
  Yvec = data[,as.character(attributes(terms(formula))$variables)[2]]
  Xmat = model.matrix(formula,data=data)
  reXmat = Xmat
  colnames(reXmat) = c('(Intercept)', 'x') 
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
  

  
  if (all(nugget == 0)){
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
  }  else {
  # $propNugget

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
  if (reml){
  thesummary$summary = summaryGmrfFit(thel)$reml
  }else{
    thesummary$summary = summaryGmrfFit(thel)$ml  
  }

  thesummary$param = thesummary$summary[,'mle']
  return(thesummary)
}

)


