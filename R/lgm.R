setGeneric('lgm', function(formula,data, ...) standardGeneric("lgm"))

setMethod("lgm", 
		signature("ANY", "Spatial"), 
		
		function(formula,  data,  
		covariates=NULL, 
		shape=1, boxcox=1, nugget = 0, 
		newdata, 
		expPred=FALSE, nuggetInPrediction=TRUE,
		reml=TRUE,mc.cores=1,
		aniso=FALSE,
		fixShape=TRUE,
		fixBoxcox=TRUE,
		fixNugget = FALSE,
		...){
	
	
	locations = newdata
	
# the formula
	# get rid of special character is names of data
	names(data) = gsub("[[:punct:]]|[[:space:]]","_", names(data))
	
	if(is.null(formula))
		formula = names(data)[1]
	if(is.integer(formula))
		formula = names(data)[formula]
	if(class(formula)!= "formula") {
		if(length(covariates)) {
			if(!length(names(covariates))) {
				names(covariates) = 
						paste("c", 1:length(covariates),sep="")			
			}
			names(covariates) = 
					gsub("[[:punct:]]|[[:space:]]",
							"_", names(covariates))

			formula = as.formula(
					paste(formula, "~ ",
							paste(names(covariates),collapse="+")
					)
			)
		} else { # end covariates not null
			formula = as.formula(paste(formula, "~1"))	
		}
	} # end formula not a formula
	
# extract covariates	
	
	# check for factors
	allterms = rownames(attributes(terms(formula))$factors)
	
	allterms = gsub("^offset\\(", "", allterms)
	alltermsWithF = gsub("\\)$", "", allterms)
	theFactors = grep("^factor", alltermsWithF, value=T)
	theFactors = gsub("^factor\\(", "", theFactors)
	
	allterms = gsub("^factor\\(", "", alltermsWithF)
	
	
	notInData = allterms[! allterms %in% names(data)]

	# convert covariates to raster stack with same resolution of prediction raster.

	if(length(covariates)){
		# extract covariate values and put in the dataset
		
		if(length(notInData)==1 & length(names(covariates))==1){
			names(covariates) = notInData
		}
		
		for(D in notInData) {
			
			if(!.compareCRS(covariates[[D]], data,unknown=TRUE) ) {
				
				require('rgdal', quietly=TRUE ) 
				
				data[[D]] = raster::extract(covariates[[D]], 
					spTransform(data, CRSobj=CRS(projection(covariates[[D]])))) 
			} else {
				data[[D]] = raster::extract(covariates[[D]], 
						 data) 
			}
			# check for factors
			
			if(!is.null(levels(covariates[[D]]))){
				# create factor, make most common value the baseline
				theTable = sort(table(data[[D]]), decreasing=TRUE)
				data[[D]] = factor(data[[D]], levels=as.integer(names(theTable)))
			}
		}	
		
		

	} 
	
	if(! all(notInData %in% names(covariates)))
		warning("some terms in the model are missing from both the data and the covariates")

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


	
# call krige	
	krigeRes =  krigeLgm(data=data,formula=formula,
			param=likRes$param, newdata=locations,
			covariates=covariates, expPred=expPred,
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