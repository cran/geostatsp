setGeneric('lgm', function(
        formula,data,grid=NULL,covariates=NULL, ...) 
			standardGeneric("lgm")
)

# sort out formula
# null formula
setMethod("lgm", 
		signature("NULL"), 
    function(formula=NULL, data, grid, 
        covariates=NULL, ...) {
    formula =  1 
  callGeneric(formula, data, grid, covariates, ...)
  }
  )


setMethod("lgm", 
		signature("numeric"),
    function(formula, data, grid, 
        covariates=NULL, ...) {
      
      formula = names(data)[formula]
      callGeneric(formula, data, grid, covariates, ...)
    }
)

# change character to formula
setMethod("lgm", 
		signature("character"),  
		function(formula, data, grid, 
        covariates=NULL, ...) {
      
      if(length(names(covariates)))
        names(covariates) = gsub("[[:punct:]]|[[:space:]]","_", names(covariates))
      if(length(covariates) & !length(names(covariates))) 
        names(covariates) = paste("c", 1:length(covariates),sep="")			
      
      if(length(formula)==1)
        formula = unique(c(formula, names(covariates)))
      if(length(formula)==1)
        formula = c(formula, '1')
      
      formula = paste(formula[1] , "~",
          paste(formula[-1], collapse=" + ")
      )
      formula = as.formula(formula)
      
      callGeneric(formula, data, grid, covariates, ...)
    }
)


# missing covariates, create empty list
setMethod("lgm", 
    signature("formula", "Spatial", "ANY", "missing"),
    function(formula, data, grid=NULL, covariates=NULL, ...) {
      callGeneric(formula, data, grid, 
          covariates=list(), 
          ...)
    }
)

# data is a raster.  grid is ignored
setMethod("lgm", 
    signature("formula", "Raster", "ANY", "ANY"),
    function(
        formula, 
        data,  
        grid=NULL,
        covariates=NULL, ...) {
      
      dataCov = gm.dataRaster(
          formula, data,
          grid=raster(data),
          covariates=covariates,
          buffer=0)
      
      callGeneric(formula, 
          dataCov$data, dataCov$grid, 
          dataCov$covariates, ...)
    }
)


# numeric cells, create raster from data bounding box

setMethod("lgm", 
		signature("formula", "Spatial", "numeric", "ANY"),
    function(formula, data, grid, covariates=NULL, ...) {
      grid = squareRaster(data, grid)
      callGeneric(formula, data, grid, covariates, ...)
    }
)


setMethod("lgm",
			signature("formula", "Spatial", "Raster", "ANY"),
      function(formula, 
          data, grid, 
          covariates=NULL, 
          buffer=0,...) {

        dataCov = gm.dataSpatial(
            formula, data, 
            grid, covariates, buffer)

        callGeneric(formula, 
            dataCov$data, dataCov$grid, 
            dataCov$covariates, ...)
      }
)

# the real work
setMethod("lgm", 
		signature("formula", "Spatial", "Raster","data.frame"), 
		function(formula, 
				data,
				grid,
		covariates=NULL, 
		shape=1, boxcox=1, nugget = 0, 
		expPred=FALSE, nuggetInPrediction=TRUE,
		reml=TRUE,mc.cores=1,
		aniso=FALSE,
		fixShape=TRUE,
		fixBoxcox=TRUE,
		fixNugget = FALSE,
		...){
  
	locations = grid
	
	
	dots <- list(...)  
	if(any(names(dots)=='param')) {
		param=dots$param	
	} else {
		param=c()
	}
	

	paramToEstimate	= c("range", "shape","nugget","boxcox")[
			!c(FALSE,fixShape,fixNugget,fixBoxcox)]		
	range=NA
	Spar = c(shape=shape,nugget=nugget,range=NA,boxcox=boxcox)
	
	if(aniso) {
		Spar = c(Spar, anisoAngleDegrees=NA,anisoRatio=NA)
		paramToEstimate = c(paramToEstimate,
				"anisoAngleDegrees","anisoRatio")		
	}

	Spar = Spar[!names(Spar) %in% names(param)]
	param = c(param, Spar)
	
	
	
	
	# to do: make sure factors in rasters are set up correctly
	# have baseline as first entry in cov@data@attributes,
	# NA's for levels without data
	# have most common level the baseline
	
# call likfit
	
	dots$param = param
	dots$formula=formula
	dots$data=data
	dots$paramToEstimate=paramToEstimate

 	likRes = do.call(likfitLgm, dots)
 
#	stuff <<- list(formula=formula, data=data, grid=grid, covariates=covariates,
#			param=likRes$param,expPred=expPred,nuggetInPrediction=nuggetInPrediction)	


# call krige	
	krigeRes =  krigeLgm(
			formula=formula,data=data,
			grid=grid,
			covariates=covariates, param=likRes$param, 
			expPred=expPred,
			nuggetInPrediction=nuggetInPrediction
			)
			
		#	data$resid = likRes$resid$resid
#	likRes$data = data

	data@data = cbind(data.frame(data), 
			likRes$data[rownames(data.frame(data)),])
	likRes$data = data
		

	res = c(predict=krigeRes, likRes)

	
	# add confidence intervals for covariance parameters
	theInf=informationLgm(res)
	res$varBetaHat = list(beta=res$varBetaHat)
	names(res) = gsub("varBetaHat", "varParam", names(res))
	res$summary = 	theInf$summary
	res$varParam$information = theInf$information

	if(FALSE){
	for(Dvar in names(covariates)) {
		theLevels =levels(covariates[[Dvar]])[[1]]
		if(!is.null(nrow(theLevels))){
			for(D in 1:nrow(theLevels)) {
				rownames(res$summary) = gsub(
						paste("(factor)?(\\()?", Dvar, "(\\))?:?", 
								theLevels[D,1],"$",sep=""),
						paste(Dvar, ":",theLevels[D,2],sep=""), 
						rownames(res$summary))
			}
		}
	}
	}
	# if range is very big, it's probably in metres, convert to km
	if(res$summary['range','estimate']>1000) {
		logicalCol = names(res$summary) == "Estimated"
		res$summary["range",!logicalCol] = 
				res$summary["range",!logicalCol] /1000
		rownames(res$summary) = gsub("^range$", "range/1000", 
				rownames(res$summary))
	}
	
	return(res)
}

)