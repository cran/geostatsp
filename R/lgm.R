lgm <- function(data,  cells, covariates=NULL, formula=NULL,
		rough=1, fixRough=TRUE,
		aniso=FALSE, boxcox=1, fixBoxcox=TRUE,
		nugget = 0, fixNugget = FALSE, ...){

# create raster
	if(!length(grep("^Raster",class(cells)))) { 
		# cells must be an integer
		cells = as.integer(cells)
		thebbox = data@bbox
		res = diff(thebbox[1,])/cells		
		Nx = cells
		Ny = ceiling(diff(thebbox[2,])/res)
		thebbox[2,2] = thebbox[2,1] + res*Ny
		cells= raster(extent(thebbox), ncols=Nx, nrows=Ny, 
				crs=data@proj4string)
	} 
	
	if(cells@nrows * cells@ncols > 10^6) warning("there are lots of cells in the prediction raster,\n this might take a very long time")
	
	
	
# the formula
	# get rid of special character is names of data
	names(data) = gsub("[[:punct:]]|[[:space:]]","_", names(data))
	
	if(is.null(formula))
		formula = names(data)[1]
	if(is.integer(formula))
		formula = names(data)[formula]
	if(class(formula)!= "formula") {
		if(!is.null(covariates)) {
			if(!length(names(covariates)))
				names(covariates) = paste("c", 1:length(covariates),sep="")			
			names(covariates) = gsub("[[:punct:]]|[[:space:]]","_", names(covariates))
			
			formula = as.formula(
					paste(formula, "~ ",
							paste(names(covariates),collapse="+")
					)
			)
		} else {
			formula = as.formula(paste(formula, "~1"))	
		}
	}
	
# extract covariates	
	
	# check for factors
	allterms = rownames(attributes(terms(formula))$factors)
	
	allterms = gsub("^offset\\(", "", allterms)
	alltermsWithF = gsub("\\)$", "", allterms)
	theFactors = grep("^factor", alltermsWithF, value=T)
	theFactors = gsub("^factor\\(", "", theFactors)
	
	notInData = allterms[! allterms %in% names(data)]
	
	# convert covariates to raster stack with same resolution of prediction raster.

	if(!is.null(covariates)){
		# extract covariate values and put in the dataset
		
		for(D in notInData) {
			data[[D]] = extract(covariates[[D]], data)
			# check for factors
			
			if(!is.null(levels(covariates[[D]]))){
				data[[D]] = factor(data[[D]])
			}
		}	
		
		
		method = rep("bilinear", length(covariates))
		
		names(method)=names(covariates)
		
		if(!all(theFactors %in% c(names(covariates), names(data)))) {
			warning("some covariates in the model aren't in the data")
		}
		method[names(method) %in% theFactors] = "ngb" 
		
		covariates = stackRasterList(covariates, cells, method=method)
		
	} 
	
	if(! all(notInData %in% names(covariates)))
		warning("some terms in the model are missing from both the data and the covariates")
	
	dotdotnames = names(list(...))
	# check for kappa and rough both specified
	if(length(grep("^kappa$", dotdotnames)) )
		warning("specify `rough' instead of `kappa'")
	if(length(grep("^fix.kappa$", dotdotnames)) )
		warning("specify `fixRough' instead of `fix.kappa'")
	if(length(grep("^cov.model$", dotdotnames)) )
		warning("do not specify cov.model, `matern' will be used.")
	if(length(grep("^(fix.psiA|fix.psiR)$", dotdotnames)) )
		warning("do not psiA or psiR, use aniso=TRUE ")
	
	
	# call likfit
	
	likRes =  likfit(geodata=data, formula=formula, 
			cov.model="matern", kappa = rough, 
		 	fix.kappa=fixRough,
			lambda=boxcox, fix.lambda=fixBoxcox,
			fix.psiA=!aniso, fix.psiR=!aniso,
			nugget=nugget, fix.nugget=fixNugget,
			...)
	
# call krige
	
 
	krigeRes =  krige(obj.model=likRes, geodata=data, 
			locations=cells, covariates=covariates,  ...
			)
 
	 
# prepare results	
	
	if(is.list(krigeRes)) {
		res = krigeRes
	} else {
		res = list(predict=krigeRes)
	}
	
 	res$likfit = likRes
	res$parameters = as.data.frame(likRes$beta.table)
	res$parameters$fixed=FALSE
	
	fill = rep(NA, dim(res$parameters)[2]-2)
	res$parameters = rbind(
			res$parameters,
			sd=c(sqrt(likRes$sigmasq), fill, 
					FALSE)
			)

	res$parameters = rbind(
				res$parameters,
				nuggetSd=c(sqrt(likRes$nugget), fill, fixNugget)
	) 
	
	res$parameters = rbind(
				res$parameters,
				Matern=c(likRes$kappa, fill, fixRough)
	)
	
	if(aniso) {
		res$parameters = rbind(
				res$parameters,
				angle=c(likRes$aniso.pars["psiA"], fill, FALSE),
				ratio = c(likRes$aniso.pars["psiR"], fill, FALSE)
		)
	}		
	
	if(likRes$lambda != 1) {
		res$parameters = rbind(
				res$parameters,
				BoxCox=c(likRes$lambda, fill, fixBoxcox)
		)
		
	}
		

	res$parameters$fixed = as.logical(res$parameters$fixed)
	
	return(res)


	
}

