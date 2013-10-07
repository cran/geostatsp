lgm <- function(data,  locations, covariates=NULL, formula=NULL,
		rough=1, fixRough=TRUE,
		aniso=FALSE, boxcox=1, fixBoxcox=TRUE,
		nugget = 0, fixNugget = FALSE,
		expPred=FALSE, nugget.in.prediction=TRUE){
	
	
# the formula
	# get rid of special character is names of data
	names(data) = gsub("[[:punct:]]|[[:space:]]","_", names(data))
	
	if(is.null(formula))
		formula = names(data)[1]
	if(is.integer(formula))
		formula = names(data)[formula]
	if(class(formula)!= "formula") {
		if(length(covariates)) {
			if(!length(names(covariates)))
				names(covariates) = paste("c", 1:length(covariates),sep="")			
			names(covariates) = gsub("[[:punct:]]|[[:space:]]","_", names(covariates))
			
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
		
		if(length(notInData)==1){
			names(covariates) = notInData
		}
		
		for(D in notInData) {
			
			if(proj4string(covariates[[D]])!= "NA" & 
					!is.na(proj4string(data)) ) {
				data[[D]] = extract(covariates[[D]], 
					spTransform(data, CRS(proj4string(covariates[[D]])))) 
			} else {
				data[[D]] = extract(covariates[[D]], 
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

	
	param = c(range=sd(data@coords[,1]),
				rough=rough, nugget=nugget,boxcox=boxcox
				)
	paramToEstimate	= c("range", "rough","nugget","boxcox")[
			!c(FALSE,fixRough,fixNugget,fixBoxcox)]		
	if(aniso) {
		param = c(param, aniso.angle.degrees=0,aniso.ratio=1)
		paramToEstimate = c(paramToEstimate,c("aniso.angle.degrees","aniso.ratio"))		
	}
				
	
	# to do: make sure factors in rasters are set up correctly
	# have baseline as first entry in cov@data@attributes,
	# NA's for levels without data
	# have most common level the baseline
	
# call likfit
	
 	
	likRes = likfitLgm(data=data, trend=formula,
			param=param, paramToEstimate=paramToEstimate
	)
	
# call krige	
	
 
	krigeRes =  krige(data=data,trend=formula,
			param=likRes$param, locations=locations,
			covariates=covariates, expPred=expPred,
			nugget.in.prediction=nugget.in.prediction
			)
	 
	res = c(predict=krigeRes, likRes)
	
	
	return(res)


	
}

