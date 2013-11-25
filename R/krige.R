krige = function(data, trend, 
		coordinates=data,
		param,  locations, covariates=NULL, 
		expPred=FALSE,
		nuggetInPrediction=TRUE) {
	
	haveBoxCox = any(names(param)=="boxcox")
	if(haveBoxCox)
		haveBoxCox = abs(param["boxcox"]-1) > 0.001
	NsimBoxCox=40
	
	if(is.numeric(locations)){
		# locations is number of cells in the x direction
		Nx = locations
		Ny = round(locations*diff(data@bbox[2,])/diff(data@bbox[1,]))
		myExtent = 	extent(data@bbox)
		myExtent@ymax = myExtent@ymin + Ny * diff(data@bbox[1,])/Nx
		locations = raster(myExtent, Ny, Nx,
				data@proj4string)	
	}
	if(nrow(locations) * ncol(locations) > 10^7) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
	
	
	
	# find factors, so we reproject rasters using
	# the correct method.
	# search for factors in the data supplied
	
	# look for factors in the model formula
	if(class(trend)=="formula"){
		trendFormula = trend		
		covariatesForData = as.data.frame(data)
		
		if(is.vector(data)) {
			observations = data
		} else {
			
			observations = formulaLhs(trend)
			
			observations = covariatesForData[,observations]
		}
		
		theVars = all.vars(formulaRhs(trendFormula))
		
		if(length(theVars)) {
			factorsInData = unlist(lapply(
							covariatesForData[,theVars,drop=FALSE],
							is.factor))
			factorsInData = names(factorsInData)[factorsInData]
		} else {
			factorsInData=NULL
		}
		
		
		allterms = rownames(attributes(terms(trend))$factors)
		
		factorsInFormula = grep("^factor", allterms, value=TRUE)
		factorsInFormula = gsub("^factor\\(", "", factorsInFormula)
		factorsInFormula = gsub("\\)$", "", factorsInFormula)
		
		factorsInTrend=NULL
		trendFormula = trend
		
		allterms = gsub("^[[:alnum:]]+\\(", "", allterms)
		allterms = gsub("\\)$", "", allterms)
		
		if(!all(allterms %in% names(data)))
			warning("some covariates don't appear in data")
	} else {
		# trend is a data frame of covariates
		# look for factors in it
		covariatesForData = as.data.frame(trend)
		
		observations = as.data.frame(data)[,1]
		
		factorsInTrend = unlist(lapply(
						covariatesForData, is.factor
				))
		factorsInTrend = names(factorsInTrend)[factorsInTrend]
		factorsInFormula = factorsInData = NULL
		
		# guess at the formula
		trendFormlua = as.formula(paste("~",
						paste(names(covariatesForData), collapse="+")
				))
		
	} # end trend not a formula
	
	
	# we know which variables factors
	theFactors = unique(c(factorsInFormula, factorsInData, factorsInTrend))
	theFactors = theFactors[theFactors %in% names(covariates) ]
	
	
	# if there's only variable in the model assign it's name to covariates
	covariateNames = all.vars(formulaRhs(trendFormula))
	if(length(covariateNames)==1){
		# so far only one variable
		names(covariates)= covariateNames
	} 
	
	
	# loop through factors
	# and make sure integer values in rasters get converted
	# to things with parameter values!
	for(D in theFactors) {
		# is this variable in param with  a factor around it? 
		# for instance factor(x)1 and factor(x)2 ?
		paramWithFactor = grep(
				paste("factor\\(", D, "\\)[[:digit:]]+$", sep=""),
				names(param), value=TRUE)
		paramStartWithD = grep(
				paste("^", D, ".+$", sep=""),
				names(param), value=TRUE)
		paramFactorCharacter = grep(
				paste("factor\\(", D, "\\).+$", sep=""),
				names(param), value=TRUE)
		if(length(paramWithFactor)) {
			# formula will convert to factor, don't 
			# create factor beforehand
			covariates[[D]]@data@isfactor = FALSE
			theLevels = gsub(
					paste("^factor\\(",D,"\\)",sep=""),
					"",paramWithFactor)
			theLevels = as.integer(theLevels)
			haveParams = values(covariates[[D]]) %in% theLevels
			# any value we don't have a parameter for
			# gets assigned to the baseline
			# give it -999 so it's made baseline by factor() 
			thevalues = values(covariates[[D]])
			thevalues[!haveParams] = -999
			values(covariates[[D]]) = thevalues
		} else if( length(paramStartWithD) ) {
			# not a bunch of digits, 
			# stuff like xTrees and xGrassland for covariate x and levels Trees and Grassland
			# see if these line up with 
			theLevels = gsub(paste("^", D, sep=""),"",paramStartWithD)
			levelsTable = covariates[[D]]@data@attributes[[1]]
			levelsInTable = levelsTable[,2]%in% theLevels
			if(mean(theLevels %in% levelsTable[,2]) < 0.4)
				warning("many levels appear missing in covariate", D)
			valuesInParams = levelsTable[levelsInTable,1]
			# assign baseline to values with no params
			# give them -999 with label of "0" (zero)
			haveParams = values(covariates[[D]]) %in% valuesInParams
			thevalues = values(covariates[[D]])
			thevalues[!haveParams] = -999
			values(covariates[[D]]) = thevalues
			
			levelsTable = 
					levelsTable[c(1, 1:nrow(levelsTable)),]
			levelsTable[1,1]= -999
			levelsTable[1,2] = "0"
			colnames(levelsTable)[2] = ""
			covariates[[D]]@data@attributes[[1]] =  levelsTable			
			covariates[[D]]@data@isfactor = TRUE
			
		} else if (length(paramFactorCharacter)) {
			# stuff like factor(x)Trees and factor(x)Grassland for covariate x and levels Trees and Grassland
			theLevels = gsub(paste("^factor\\(", D,"\\)", sep=""),"",
					paramFactorCharacter)
			levelsTable = covariates[[D]]@data@attributes[[1]]
			levelsInTable = levelsTable[,2]%in% theLevels
			if(mean(theLevels %in% levelsTable[,2]) < 0.4)
				warning("many levels appear missing in covariate", D)
			valuesInParams = as.numeric(levelsTable[levelsInTable,1])
			# assign baseline to values with no params
			# give them -999 with label of "0" (zero)
			thevalues = values(covariates[[D]])
			haveParams = thevalues %in% valuesInParams
			thevalues[!haveParams] = -999
			values(covariates[[D]]) = thevalues
			
			levelsTable = 
					levelsTable[c(1, 1:nrow(levelsTable)),]
			levelsTable[1,1]= -999
			levelsTable[1,2] = "0"
			colnames(levelsTable)[2]=""
			covariates[[D]]@data@attributes[[1]] =  levelsTable			
			
			
		} else {
			warning("don't know what to do with covariate", D, 
					"\n can't assign parameters to levels of this factor")			
		}
		
		
	}
	
	
	if(!is.null(covariates) & length(theVars)) {
		# method for resampling covariate rasters
		
		method = rep("bilinear", length(names(covariates)))
		names(method) = names(covariates)
		method[theFactors] = "ngb"
		
		covariates = stackRasterList(covariates, locations, method=method)
		
		theVars = all.vars(formulaRhs(trendFormula))
		if(nlayers(covariates)==1 & length(theVars)==1) {
			names(covariates) = theVars
		}
		
		# construct the fixed effects component
		covariatesDF = as.data.frame(covariates, xy=TRUE)
		# convert to factors
		# get rid of trailing _ created by as.data.frame
		names(covariatesDF) = gsub("_$", "", names(covariatesDF))
	} else {
		covariatesDF = as.data.frame(matrix(NA, ncol=0, nrow=ncell(locations)))
	}
	
	
	
	# convert trend formula to LHS
	trendFormula = formulaRhs(trendFormula)
	meanRaster = raster(locations)
	
	
	missingVars = all.vars(trendFormula)%in% names(covariatesDF)
	missingVars = all.vars(trendFormula)[!missingVars]
	
	# check if all variables are in covariates
	if(length(missingVars)) {
		cat("cant find covariates ",
				paste(missingVars, collapse=","),
				" for prediction, imputing zeros\n")		
		
		covariatesDF[,missingVars]=0	
	}
	modelMatrixForRaster = model.matrix(trendFormula, covariatesDF)
	
	meanFixedEffects = 
			modelMatrixForRaster %*% param[colnames(modelMatrixForRaster)]
	
	anyNA = apply(covariatesDF, 1, function(qq) any(is.na(qq)))
	if(any(anyNA)) {
		oldmm = rep(NA, ncell(meanRaster))
		oldmm[!anyNA] = meanFixedEffects
		values(meanRaster) = oldmm
	} else {
		values(meanRaster) = meanFixedEffects
	}
	names(meanRaster) = "fixed"
	
	
	if(!all(colnames(modelMatrixForRaster)%in% names(param))){
		warning("cant find coefficients",
				paste(names(modelMatrixForRaster)[
								!names(modelMatrixForRaster)%in% names(param)
						], collapse=","),
				"in param\n")
	}
	
	
	
	
# subtract mean from data
	theNAdata = apply(covariatesForData, 1, function(qq) any(is.na(qq))) | 
			is.na(observations)
	
	
	if(any(theNAdata)) {
		noNAdata = !theNAdata
		if(length(grep("^SpatialPoints", class(coordinates)))) {
			coordinates = SpatialPoints(coordinates)[noNAdata]	
		} else if(class(coordinates)=="dist"){
			coordinates = as.matrix(coordinates)
			coordinates = coordinates[noNAdata,noNAdata]
			coordinates = as.dist(coordinates)
		} else {
			warning("missing vlaues in data but unclear how to remove them from coordinates")
		}
		observations = observations[noNAdata]
		covariatesForData = covariatesForData[noNAdata,,drop=FALSE]
	}
	
	
	modelMatrixForData = model.matrix(trendFormula, covariatesForData)
	
	meanForData = 
			modelMatrixForData %*% param[colnames(modelMatrixForData)]
	
	
	if(haveBoxCox) {
		if(abs(param["boxcox"]) < 0.001) {
			observations = log(observations)
			expPred = TRUE
			haveBoxCox = FALSE
		} else {
			observations = ((observations^param["boxcox"]) - 1)/
					param["boxcox"]
		}
		
	}
	
	observations = observations - meanForData
	
	haveNugget = any(names(param)=="nugget")
	if(haveNugget) { 
		haveNugget = param["nugget"] > 0
	} 
	if(!haveNugget) {
		nuggetInPrediction=FALSE
	}	
	
	varData = matern(coordinates, param=param)
	if(haveNugget) Matrix::diag(varData) = Matrix::diag(varData) + param["nugget"]
	
	cholVarData = Matrix::chol(varData)
	cholVarDatInvData = Matrix::solve(cholVarData, observations)


	Ny = length(observations)
	param = fillParam(param)
	
	krigeOneRow = function(Drow){

		# covariance of cells in row Drow with data points
		resC =  .C("maternArasterBpoints", 
				as.double(xmin(locations)), as.double(xres(locations)), 
						as.integer(ncol(locations)), 
				as.double(yFromRow(locations, Drow)), 
					as.double(0), as.integer(1),
				as.double(coordinates@coords[,1]), 
					as.double(coordinates@coords[,2]), 
				N=as.integer(Ny), 
				result=as.double(matrix(0, ncol(locations), length(coordinates))),
				as.double(param["range"]),
				as.double(param["shape"]),
				as.double(param["variance"]),
				as.double(param["anisoRatio"]),
				as.double(param["anisoAngleRadians"])
		) 
		covDataPred = matrix(resC$result, nrow=ncol(locations), ncol=Ny)
		
		
		cholVarDataInvCovDataPred = Matrix::solve(cholVarData, t(covDataPred))

		c( # the conditional expectation
			as.vector(Matrix::crossprod(cholVarDataInvCovDataPred, 
						cholVarDatInvData)),
		# part of the conditional variance
			apply(cholVarDataInvCovDataPred, 2, function(qq) {
					sumsq = sum(qq^2)
				})
		)
		
	}

	sums = mcmapply(krigeOneRow,1:nrow(locations))

	# row sums of cholVarDataInvCovDataPred
	forExpected = sums[1:ncol(locations),]
	# row sums of squares
	forVar = as.vector(sums[seq(ncol(locations)+1, len=ncol(locations)),])
	
	
	randomRaster = raster(meanRaster)
	names(randomRaster) = "random"
 	values(randomRaster) = as.vector(forExpected)
	
	predRaster = meanRaster + randomRaster
	names(predRaster) = "predict"
	
	
	
	if(any(forVar > param["variance"])){
		themax = max(forVar - param["variance"],na.rm=TRUE)
		if(themax > 1e-6)
			warning("converted variances of ", themax, " to zero")	
		forVar = pmin(forVar, param["variance"])	
	}
	
	krigeSd = raster(meanRaster)
	names(krigeSd) = "krigeSd"

	if(nuggetInPrediction) {
		values(krigeSd) = sqrt(sum(param[c("nugget","variance")]) - 
						forVar)
	} else {
		values(krigeSd) = sqrt(param["variance"] - forVar)
	}
	
	result = stack(meanRaster, randomRaster, predRaster,
			krigeSd)

	# box-cox
	if(haveBoxCox){ 
		names(result)[names(result)=="predict"] = "predict.boxcox"
		
		themean = values(result[["predict.boxcox"]])
		thesd = values(result[["krigeSd"]])
		
		theNA = is.na(themean)
		thesd[is.na(thesd)] = 0 
		themean[is.na(themean)] = 0 
		
		invlambda = 1/param["boxcox"]
		
		
		Ndata = length(themean)
		
		# if lambda is fractional, truncate transformed values at zero
		if(is.nan((-1)^param["boxcox"])) {
			useMax=0
		} else {
			useMax = -Inf
		}
		
		bcpred = 0
		for(D in 1:NsimBoxCox) {
			bcpred = bcpred + 
					exp(invlambda*log(
									pmax(param["boxcox"] *rnorm(Ndata, themean, thesd)+1,useMax)))
		}
		bcpred[is.na(themean)] = NA
		bcpred = bcpred / NsimBoxCox
		bcpred[theNA] = NA
		
		newraster=raster(result[["predict.boxcox"]])
		names(newraster) = "predict"
		values(newraster) = bcpred
		
		result = addLayer(result, 
				newraster)
		
	} # end have box cox
	
	
	if(expPred){
		
		names(result)[names(result)=="predict"] = "predict.log"
		newLayer = exp(result[["predict.log"]]+ 0.5*result[["krigeSd"]]^2 )
		names(newLayer) = "predict"
		result = addLayer(result, newLayer)
		
	} # end expPred
		
	result
}

krigeOld = function(data, trend, 
		coordinates=data,
		param,  locations, covariates=NULL, 
		expPred=FALSE,
		nuggetInPrediction=TRUE) {

	haveBoxCox = any(names(param)=="boxcox")
	if(haveBoxCox)
		haveBoxCox = abs(param["boxcox"]-1) > 0.001
	NsimBoxCox=40

	if(is.numeric(locations)){
		# locations is number of cells in the x direction
		Nx = locations
		Ny = round(locations*diff(data@bbox[2,])/diff(data@bbox[1,]))
		myExtent = 	extent(data@bbox)
		myExtent@ymax = myExtent@ymin + Ny * diff(data@bbox[1,])/Nx
		locations = raster(myExtent, Ny, Nx,
			 data@proj4string)	
	}
	if(nrow(locations) * ncol(locations) > 10^7) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
		
		
	
	# find factors, so we reproject rasters using
	# the correct method.
	# search for factors in the data supplied
	
	# look for factors in the model formula
	if(class(trend)=="formula"){
		trendFormula = trend		
		covariatesForData = as.data.frame(data)
		
		if(is.vector(data)) {
			observations = data
		} else {

			observations = formulaLhs(trend)

			observations = covariatesForData[,observations]
		}
		
		theVars = all.vars(formulaRhs(trendFormula))
		
		if(length(theVars)) {
			factorsInData = unlist(lapply(
						covariatesForData[,theVars,drop=FALSE],
						is.factor))
			factorsInData = names(factorsInData)[factorsInData]
		} else {
			factorsInData=NULL
		}
		
		
		allterms = rownames(attributes(terms(trend))$factors)
	
		factorsInFormula = grep("^factor", allterms, value=TRUE)
		factorsInFormula = gsub("^factor\\(", "", factorsInFormula)
		factorsInFormula = gsub("\\)$", "", factorsInFormula)

		factorsInTrend=NULL
		trendFormula = trend
		
		allterms = gsub("^[[:alnum:]]+\\(", "", allterms)
		allterms = gsub("\\)$", "", allterms)
		
		if(!all(allterms %in% names(data)))
			warning("some covariates don't appear in data")
	} else {
		# trend is a data frame of covariates
		# look for factors in it
		covariatesForData = as.data.frame(trend)
		
		observations = as.data.frame(data)[,1]
		
		factorsInTrend = unlist(lapply(
						covariatesForData, is.factor
						))
		factorsInTrend = names(factorsInTrend)[factorsInTrend]
		factorsInFormula = factorsInData = NULL
				
		# guess at the formula
		trendFormlua = as.formula(paste("~",
						paste(names(covariatesForData), collapse="+")
						))
		
	} # end trend not a formula

		
	# we know which variables factors
	theFactors = unique(c(factorsInFormula, factorsInData, factorsInTrend))
	theFactors = theFactors[theFactors %in% names(covariates) ]
	
	
	# if there's only variable in the model assign it's name to covariates
	covariateNames = all.vars(formulaRhs(trendFormula))
	if(length(covariateNames)==1){
		# so far only one variable
		names(covariates)= covariateNames
	} 
	
	
	# loop through factors
	# and make sure integer values in rasters get converted
	# to things with parameter values!
	for(D in theFactors) {
		# is this variable in param with  a factor around it? 
		# for instance factor(x)1 and factor(x)2 ?
		paramWithFactor = grep(
				paste("factor\\(", D, "\\)[[:digit:]]+$", sep=""),
				names(param), value=TRUE)
		paramStartWithD = grep(
				paste("^", D, ".+$", sep=""),
				names(param), value=TRUE)
		paramFactorCharacter = grep(
				paste("factor\\(", D, "\\).+$", sep=""),
				names(param), value=TRUE)
 		if(length(paramWithFactor)) {
			# formula will convert to factor, don't 
			# create factor beforehand
			covariates[[D]]@data@isfactor = FALSE
			theLevels = gsub(
					paste("^factor\\(",D,"\\)",sep=""),
					"",paramWithFactor)
			theLevels = as.integer(theLevels)
			haveParams = values(covariates[[D]]) %in% theLevels
			# any value we don't have a parameter for
			# gets assigned to the baseline
			# give it -999 so it's made baseline by factor() 
			thevalues = values(covariates[[D]])
			thevalues[!haveParams] = -999
			values(covariates[[D]]) = thevalues
		} else if( length(paramStartWithD) ) {
			# not a bunch of digits, 
			# stuff like xTrees and xGrassland for covariate x and levels Trees and Grassland
			# see if these line up with 
			theLevels = gsub(paste("^", D, sep=""),"",paramStartWithD)
			levelsTable = covariates[[D]]@data@attributes[[1]]
			levelsInTable = levelsTable[,2]%in% theLevels
			if(mean(theLevels %in% levelsTable[,2]) < 0.4)
				warning("many levels appear missing in covariate", D)
			valuesInParams = levelsTable[levelsInTable,1]
			# assign baseline to values with no params
			# give them -999 with label of "0" (zero)
			haveParams = values(covariates[[D]]) %in% valuesInParams
			thevalues = values(covariates[[D]])
			thevalues[!haveParams] = -999
			values(covariates[[D]]) = thevalues
		
			levelsTable = 
					levelsTable[c(1, 1:nrow(levelsTable)),]
			levelsTable[1,1]= -999
			levelsTable[1,2] = "0"
			colnames(levelsTable)[2] = ""
			covariates[[D]]@data@attributes[[1]] =  levelsTable			
			covariates[[D]]@data@isfactor = TRUE
			
		} else if (length(paramFactorCharacter)) {
			# stuff like factor(x)Trees and factor(x)Grassland for covariate x and levels Trees and Grassland
	theLevels = gsub(paste("^factor\\(", D,"\\)", sep=""),"",
			paramFactorCharacter)
	levelsTable = covariates[[D]]@data@attributes[[1]]
	levelsInTable = levelsTable[,2]%in% theLevels
	if(mean(theLevels %in% levelsTable[,2]) < 0.4)
		warning("many levels appear missing in covariate", D)
	valuesInParams = as.numeric(levelsTable[levelsInTable,1])
	# assign baseline to values with no params
	# give them -999 with label of "0" (zero)
	thevalues = values(covariates[[D]])
	haveParams = thevalues %in% valuesInParams
	thevalues[!haveParams] = -999
	values(covariates[[D]]) = thevalues
	
	levelsTable = 
			levelsTable[c(1, 1:nrow(levelsTable)),]
	levelsTable[1,1]= -999
	levelsTable[1,2] = "0"
	colnames(levelsTable)[2]=""
	covariates[[D]]@data@attributes[[1]] =  levelsTable			
	
			
		} else {
			warning("don't know what to do with covariate", D, 
					"\n can't assign parameters to levels of this factor")			
		}
		
		
	}
	
	
	if(!is.null(covariates) & length(theVars)) {
		# method for resampling covariate rasters

		method = rep("bilinear", length(names(covariates)))
		names(method) = names(covariates)
		method[theFactors] = "ngb"
		
		covariates = stackRasterList(covariates, locations, method=method)

		theVars = all.vars(formulaRhs(trendFormula))
		if(nlayers(covariates)==1 & length(theVars)==1) {
			names(covariates) = theVars
		}
		
		# construct the fixed effects component
		covariatesDF = as.data.frame(covariates, xy=TRUE)
		# convert to factors
		# get rid of trailing _ created by as.data.frame
		names(covariatesDF) = gsub("_$", "", names(covariatesDF))
	} else {
		covariatesDF = as.data.frame(matrix(NA, ncol=0, nrow=ncell(locations)))
	}
		

	
	# convert trend formula to LHS
	trendFormula = formulaRhs(trendFormula)
	meanRaster = raster(locations)

	
	missingVars = all.vars(trendFormula)%in% names(covariatesDF)
	missingVars = all.vars(trendFormula)[!missingVars]
	
	# check if all variables are in covariates
	if(length(missingVars)) {
		cat("cant find covariates ",
				paste(missingVars, collapse=","),
				" for prediction, imputing zeros\n")		
		
		covariatesDF[,missingVars]=0	
	}
	modelMatrixForRaster = model.matrix(trendFormula, covariatesDF)

	meanFixedEffects = 
			modelMatrixForRaster %*% param[colnames(modelMatrixForRaster)]

	anyNA = apply(covariatesDF, 1, function(qq) any(is.na(qq)))
	if(any(anyNA)) {
		oldmm = rep(NA, ncell(meanRaster))
		oldmm[!anyNA] = meanFixedEffects
		values(meanRaster) = oldmm
	} else {
		values(meanRaster) = meanFixedEffects
	}
	names(meanRaster) = "fixed"


	if(!all(colnames(modelMatrixForRaster)%in% names(param))){
		warning("cant find coefficients",
				paste(names(modelMatrixForRaster)[
								!names(modelMatrixForRaster)%in% names(param)
						], collapse=","),
				"in param\n")
	}
	
	
	
	
# subtract mean from data
theNAdata = apply(covariatesForData, 1, function(qq) any(is.na(qq))) | 
		is.na(observations)


if(any(theNAdata)) {
	noNAdata = !theNAdata
	if(length(grep("^SpatialPoints", class(coordinates)))) {
		coordinates = SpatialPoints(coordinates)[noNAdata]	
	} else if(class(coordinates)=="dist"){
		coordinates = as.matrix(coordinates)
		coordinates = coordinates[noNAdata,noNAdata]
		coordinates = as.dist(coordinates)
	} else {
		warning("missing vlaues in data but unclear how to remove them from coordinates")
	}
	observations = observations[noNAdata]
	covariatesForData = covariatesForData[noNAdata,,drop=FALSE]
}


modelMatrixForData = model.matrix(trendFormula, covariatesForData)

meanForData = 
		modelMatrixForData %*% param[colnames(modelMatrixForData)]


if(haveBoxCox) {
	if(abs(param["boxcox"]) < 0.001) {
		observations = log(observations)
		expPred = TRUE
		haveBoxCox = FALSE
	} else {
		observations = ((observations^param["boxcox"]) - 1)/
				param["boxcox"]
	}
	
}

observations = observations - meanForData

haveNugget = any(names(param)=="nugget")
if(haveNugget) { 
	haveNugget = param["nugget"] > 0
} 
if(!haveNugget) {
	nuggetInPrediction=FALSE
}	
	
varData = matern(coordinates, param=param)
if(haveNugget) Matrix::diag(varData) = Matrix::diag(varData) + param["nugget"]

cholVarData = Matrix::chol(varData)
cholVarDatInvData = Matrix::solve(cholVarData, observations)

covDataPred = matern(coordinates,locations, param=param)
cholVarDataInvCovDataPred = Matrix::solve(cholVarData, covDataPred)

randomRaster = raster(meanRaster)
values(randomRaster) = Matrix::crossprod(cholVarDataInvCovDataPred, cholVarDatInvData)[,1]
names(randomRaster) = "random"

predRaster = meanRaster + randomRaster
names(predRaster) = "predict"


krigeSd = raster(meanRaster)
thesd = apply(cholVarDataInvCovDataPred, 2,function(qq) sum(qq*qq))

# check for values bigger than the variance
if(any(thesd > param["variance"])){
	themax = max(thesd - param["variance"])
	if(themax > 1e-6)
		warning("converted variances of ", themax, " to zero")	
	thesd = pmin(thesd, param["variance"])	
}

values(krigeSd) = thesd
 # krigeSd is now a variance
# transform to Standard Deviation

if(nuggetInPrediction) {
	values(krigeSd) = sqrt(sum(param[c("nugget","variance")]) - 
					values(krigeSd))
} else {
	values(krigeSd) = sqrt(param["variance"] - values(krigeSd))
}
	
names(krigeSd) = "krigeSd"


result = stack(meanRaster, randomRaster, predRaster,
		krigeSd)




# box-cox
if(haveBoxCox){ 
	 
		
	names(result)[names(result)=="predict"] = "predict.boxcox"

	themean = values(result[["predict.boxcox"]])
	thesd = values(result[["krigeSd"]])

	theNA = is.na(themean)
	thesd[is.na(thesd)] = 0 
	themean[is.na(themean)] = 0 
			
	invlambda = 1/param["boxcox"]
	

	Ndata = length(themean)
	
	# if lambda is fractional, truncate transformed values at zero
	if(is.nan((-1)^param["boxcox"])) {
		useMax=0
	} else {
		useMax = -Inf
	}
	
	bcpred = 0
	for(D in 1:NsimBoxCox) {
		bcpred = bcpred + 
				exp(invlambda*log(
		pmax(param["boxcox"] *rnorm(Ndata, themean, thesd)+1,useMax)))
	}
	bcpred[is.na(themean)] = NA
	bcpred = bcpred / NsimBoxCox
	bcpred[theNA] = NA
	
	newraster=raster(result[["predict.boxcox"]])
	names(newraster) = "predict"
	values(newraster) = bcpred
	
	result = addLayer(result, 
			newraster)
	
} # end have box cox


if(expPred){
	
	names(result)[names(result)=="predict"] = "predict.log"
	newLayer = exp(result[["predict.log"]]+ 0.5*result[["krigeSd"]]^2 )
	names(newLayer) = "predict"
	result = addLayer(result, newLayer)
	
} # end expPred

conditionalVarianceMatrix = FALSE
if(conditionalVarianceMatrix) {
	condVarFull = matern(result,  param=param)
	condVar = condVarFull - Matrix::crossprod(cholVarDataInvCovDataPred)
	
	result = list(raster=result,
			condVar =  condVar)#, cholmat = cholVarDataInvCovDataPred,
#			full=condVarFull)
	
}

result
}