simPoissonPP = function(intensity) {
	
	NperCell = intensity
	values(NperCell) = rpois(ncell(intensity), 
			values(intensity)*prod(res(intensity)))
	
	events = rep(1:ncell(NperCell), values(NperCell))
	
	events = as.data.frame(NperCell,xy=TRUE)[events, c("x","y")]
	
	events = events + cbind(
			runif(dim(events)[1],-xres(intensity)/2, xres(intensity)/2),
			runif(dim(events)[1],-yres(intensity)/2, yres(intensity)/2)
	)
	
	events = SpatialPoints(events)
	
	if(!is.na(intensity@crs@projargs))
		events@proj4string = intensity@crs
	
	events
	
}

simLgcp = function(param, covariates=NULL, betas=NULL, 
		rasterTemplate=covariates[[1]], model="whittle", ...) {
	
	param["nugget"] = 0
	param = param[c("mean","variance","nugget", 
				"scale" ,"alpha")]
	
	if(!is.null(covariates))
		covariates = stackRasterList(covariates, rasterTemplate)

	randomEffect = GaussRF(rasterTemplate, model=model, 
					param=param, ...)

	linearPredictor = randomEffect
	if(is.null(names(betas)))
		names(betas) = names(covariates)
	
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