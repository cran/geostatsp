grfConditional = function(data, trend, param, locations, Nsim, 
		covariates,   fun, nugget.in.prediction=TRUE) {
	
	krigeRes = krige(data=data, trend=trend, 
			param=param, locations=locations, 
    covariates = covariates,  nugget.in.prediction = nugget.in.prediction,
	conditionalVarianceMatrix=TRUE)
	

 

locations = raster(krigeRes$raster)
Npred=ncell(locations)

condVar = krigeRes$condVar

condVarChol = Matrix::chol(condVar)

if(any(names(krigeRes$raster)=="predict.boxcox" )) {
	predPlusFit = values(krigeRes$raster[["predict.boxcox"]])
} else if(any(names(krigeRes$raster)=="predict.log") ){
	predPlusFit = values(krigeRes$raster[["predict"]])
} else {
	predPlusFit = values(krigeRes$raster[["predict"]])
}

if(missing(fun)) {
	simBig=rep(0.1, Nsim*Npred)
	for(D in 1:Nsim) {
		sim = predPlusFit  + condVarChol %*% rnorm(Npred)
		
		simBig[seq(Npred*(D-1)+1, len=Npred)] = as.vector(sim)
	}

	result = brick(array(0.1,c(nrow(locations), ncol(locations),Nsim)), 
			xmn=locations@extent@xmin, 
		xmx=locations@extent@xmax, ymn=locations@extent@ymin, 
		ymx=locations@extent@ymax,crs=locations@crs )
	values(result) = simBig	
} else { # a function was supplied
		result = list()			
		for(D in 1:Nsim) {
			result[[D]] = fun(predPlusFit + condVarChol %*% rnorm(Npred))			
		}
}

result
}



