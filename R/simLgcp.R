simPoissonPP = function(intensity) {
	
	ivec = values(intensity)
	ivec[is.na(ivec)] = 0
	ivec =  ivec*prod(res(intensity))
	
	NperCell = intensity
	values(NperCell) = rpois(ncell(intensity), ivec)
	
	if(maxValue(NperCell)>1000)
		warning("A large number of events are being simulated, more than", maxValue(NperCell))
	
	events = rep(1:ncell(NperCell), values(NperCell))
	
	if(length(events)>1e6)
		warning("more than 1,000,000 events being simulated")
	
	events = as.data.frame(NperCell,xy=TRUE)[events, c("x","y")]
	
	events = events + cbind(
			runif(dim(events)[1],-xres(intensity)/2, xres(intensity)/2),
			runif(dim(events)[1],-yres(intensity)/2, yres(intensity)/2)
	)
	
	events = SpatialPoints(events)
	projection(events) = projection(intensity)
	
	events
	
}

simLgcp = function(param, covariates=NULL, betas=NULL, 
		rasterTemplate=covariates[[1]],  ...) {
	
	randomEffect = geostatsp::RFsimulate(model=param, x=rasterTemplate,  ...)
	
	if(!is.null(covariates))
		covariates = stackRasterList(covariates, randomEffect)

	if(is.null(names(betas)))
		names(betas) = names(covariates)
	
	
	themean = 0
	if('mean' %in% names(param))
		themean = themean + param['mean']
	if('intercept' %in% names(betas))
		themean = themean + betas['intercept']
	betas['intercept'] = themean
	param = param[! names(param) %in% "mean"]

	
	
	
	linearPredictor = themean +randomEffect
	
	for(D in names(covariates)) {
		linearPredictor = linearPredictor + betas[D]* covariates[[D]]		
	}

	intensity = exp(linearPredictor)
	
	events = simPoissonPP(intensity)
	
	names(linearPredictor) = "linearPredictor"
	names(intensity) = "intensity"
	names(randomEffect) = "random"
	
	return(list(events=events,
					raster = stack(randomEffect,
							linearPredictor, intensity,covariates),
					params=list(random=param, fixed=betas)
			))
			

}