setGeneric('lgm', function(formula,data,grid,covariates=NULL, ...) standardGeneric("lgm"))

# sort out formula
# null formula
setMethod("lgm", 
		signature("NULL"), 
		gm.nullFormula
)


setMethod("lgm", 
		signature("numeric"),  
		gm.numericFormula
)

# change character to formula
setMethod("lgm", 
		signature("character"),  
		gm.characterFormula
)


# numeric cells, create raster from data bounding box

setMethod("lgm", 
		signature("formula", "ANY", "numeric"),
		gm.gridNumeric
)



# extrat covariates for data, convert covariates to a stack
setMethod("lgm", 
		signature("formula", "Raster", "Raster"),
		gm.dataRaster
)

#setMethod("lgm", 
#		signature("formula", "Raster", "missing", "ANY"),
# aggregate the covariates to data, merge covariates into data, call formula, Raster, missing, missing)
#		gm.dataRaster
#				
#)

setMethod("lgm", 
		signature("formula", "Raster", "missing", "missing"),
		function(formula, 
				data,  
				grid,
				covariates, ...) {
		gridHere = raster(data)
		if(abs(diff(res(gridHere)))>0.000001 )
			warning("data is not on a square grid")
		dataSP = as(data, "SpatialPointsDataFrame")
		dataDF = dataSP@data
		
		callGeneric(
				formula=formula, data=dataDF,
				grid=gridHere,
				covariates=data.frame(), 
				...
		)	
		}
)


setMethod("lgm", 
			signature("formula", "Spatial", "Raster", "list"),
			gm.dataSpatial
	)

	setMethod("lgm", 
			signature("formula", "Spatial", "Raster", "NULL"),
			gm.dataSpatial
	)
	
	
setMethod("lgm", 
		signature("formula", "Spatial", "Raster", "Raster"),
		gm.dataSpatial
)



setMethod("lgm", 
		signature("formula", "Spatial", "Raster","data.frame"), 
		function(formula, 
				data,  
				grid,
		covariates=list(), 
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