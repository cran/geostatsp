

krige = function(obj.model, geodata,  locations, covariates, 
		locations.mean=locations,
		factor.info=NULL, exp.pred=FALSE,rasterMethod = c("ngb", "bilinear"),
		nugget.in.prediction=TRUE, ...) {

NsimBoxCox=40

if(is.numeric(locations)){
	# locations is number of cells in the x direction
	Nx = locations
	Ny = round(locations*diff(geodata@bbox[2,])/diff(geodata@bbox[1,]))
	myExtent = 	extent(geodata@bbox)
	myExtent@ymax = myExtent@ymin + Ny * diff(geodata@bbox[1,])/Nx
	locations = raster(myExtent, Ny, Nx,
			 geodata@proj4string)	
}
if(is.numeric(locations.mean)){
	# locations is number of cells in the x direction
		Nx = locations.mean
		Ny = round(locations.mean*diff(geodata@bbox[2,])/diff(geodata@bbox[1,]))
		myExtent = 	extent(geodata@bbox)
		myExtent@ymax = myExtent@ymin + Ny * diff(geodata@bbox[1,])/Nx
		locations.mean = raster(myExtent, Ny, Nx,
				geodata@proj4string)	
	}

 
data = geodata@data	
theCovs = attributes(terms(obj.model$formula))$term.labels

varTypes = unlist(lapply(data, class))[theCovs]
thefactors = names(varTypes[varTypes=="factor"])

# if necessary, turn covariates into a raster stack
if(class(covariates) == "list") {
	
	themethod = rep(rasterMethod[1], length(covariates))
	names(themethod) = names(covariates)
	themethod[thefactors] = "ngb"
	
	locations.mean = stackRasterList(covariates,locations.mean, themethod)	
} else if (class(covariates)=="RasterStack"){ # covariates must be a raster stack, use as-is
	if(length(names(covariates))==1 & length(theCovs)==1) {
		names(covariates) = theCovs
	}
	locations.mean = covariates
} else {
	# covariates must be null
}


if(!all(names(locations.mean)%in% theCovs))
	warning("some covariates in the model formula weren't supplied\n they will be assumed to be zeros when making spatial predictions")

# data frame of random field prediction locations 
locationsDF = as.data.frame(locations, xy=TRUE)	

# convert data to factors, using same levels for model fit data and prediction data


# construct the fixed effects component
meanRaster = raster(locations.mean)
meanRaster[] = obj.model$beta["(Intercept)"]

for(D in theCovs[theCovs %in% names(covariates)]){

	if(varTypes[D] == "factor") {
		if (D %in% names(factor.info)) {
			tofac = factor.info[[D]]
			if(all(c("levels","labels") %in% names(tofac)))
				tofac=data.frame(tofac[["levels"]],tofac[["labels"]])
			rownames(tofac) = tofac[,2]
			thelevels = tofac[levels(data[[D]]),1]
		}  else {
			thelevels = levels(data[[D]])	
		}
		theBetaNames = paste(D, levels(data[[D]]), sep="")
		betaHasPar = theBetaNames %in% names(obj.model$beta)
		betasHere = rep(0, length(theBetaNames))
		names(betasHere) = theBetaNames
		betasHere[betaHasPar] = obj.model$beta[theBetaNames[betaHasPar]]
		names(betasHere) = thelevels
		
		
		toAdd = setValues(locations.mean[[D]],
				betasHere[as.character(
								locations.mean[[D]]@data@values
						)]
		)

		meanRaster = meanRaster + toAdd
		
	} else {

		meanRaster = meanRaster + obj.model$beta[D]*locations.mean[[D]]
	}
}
names(meanRaster) = "fixed"
# do the kriging

data.col = unlist(strsplit(as.character(obj.model$formula), "~"))
data.col = gsub("[[:space:]]", "", data.col)
data.col = data.col[data.col != ""]
data.col = data.col[1]
geodataForKrige = geoR::as.geodata(geodata, 
		data.col=data.col, 
		covar.col=theCovs)
if(obj.model$lambda != 1) {
	if(obj.model$lambda == 0) {
		geodataForKrige$data =log(geodataForKrige$data)
	} else{
		geodataForKrige$data =(geodataForKrige$data^obj.model$lambda -1)/obj.model$lambda
	} 
}
dummyDF = locationsDF
for(D in theCovs){
	if(varTypes[D] == "factor") {
		dummyDF[,D] = factor(levels(data[,D])[1], levels=levels(data[,D]))
	}
	else
		dummyDF[,D] = 0	
}
trend.d=  trend.spatial(obj.model$trend, data)	
trend.l =  trend.spatial(obj.model$trend, dummyDF)

obj.modelNoLambda = obj.model
obj.modelNoLambda$lambda =1
thecontrol = geoR::krige.control(obj.model=obj.modelNoLambda,
		trend.d=trend.d, 
		trend.l=trend.l) 


krigeResult = krige.conv(
		geodataForKrige, locations=	locationsDF[,1:2],
		krige=thecontrol, output = output.control(messages=FALSE)
)

rastKrige = setValues(locations[[1]], krigeResult$predict - obj.model$beta["(Intercept)"])
rastKrige = addLayer(rastKrige, 
		setValues(locations[[1]], krigeResult$krige.var)
)
names(rastKrige) = c("predict","krige.var")


if(as(meanRaster, "BasicRaster")==as(rastKrige, "BasicRaster")) {
	result = addLayer(meanRaster,
			rastKrige
	)
}	 else {
	result = addLayer(meanRaster,
			raster::resample(rastKrige, meanRaster, method="ngb")
	)
}

names(result)[names(result)=="predict"] = "random"
result = addLayer(result, 
		predict=result[["fixed"]] + result[["random"]]
)
names(result)[names(result)=="layer"] = "predict"


if(exp.pred | obj.model$lambda==0){
	names(result)[names(result)=="predict"] = "predict.log"
	result = addLayer(result, 
			exp(result[["predict.log"]]+ 0.5*result[["krige.var"]] +
							0.5* nugget.in.prediction * obj.model$nugget)
	)
	names(result)[names(result)=="layer"] = "predict"
	
}

# box-cox
if(!any(obj.model$lambda==c(0,1))){

	names(result)[names(result)=="predict"] = "predict.boxcox"

	

	
	themean = values(result[["predict.boxcox"]])

	if(nugget.in.prediction){
		thesd = sqrt(values(result[["krige.var"]]) + obj.model$nugget)
 	} else {
		thesd = sqrt(values(result[["krige.var"]]))
	}

	theNA = is.na(themean)
	thesd[is.na(thesd)] = 0 
	themean[is.na(themean)] = 0 
			
	invlambda = 1/obj.model$lambda
	
	bcpred = 0
	Ndata = length(themean)
	
	# if lambda is fractional, truncate transformed values at zero
	if(is.nan((-1)^obj.model$lambda)) {
		useMax=0
	} else {
		useMax = -Inf
	}
	
	for(D in 1:NsimBoxCox) {
		bcpred = bcpred + 
				exp(invlambda*log(
		pmax(obj.model$lambda *rnorm(Ndata, themean, thesd)+1,useMax)))
	}
	bcpred[is.na(themean)] = NA
	bcpred = bcpred / NsimBoxCox
	bcpred[theNA] = NA
	
	newraster=raster(result[["predict.boxcox"]])
	names(newraster) = "predict"
	values(newraster) = bcpred
	
	result = addLayer(result, 
			newraster)
	
}


if(as(result, 'BasicRaster')!=as(rastKrige, "BasicRaster")) {
	result = list(random = rastKrige, prediction=result)			
} 		

result
}