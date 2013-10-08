grfConditional = function(data, ycol=1, 
			param, locations, Nsim,
		 	fun=NULL, nugget.in.prediction=TRUE){
		
if(is.numeric(locations)){
	# locations is number of cells in the x direction
	Nx = locations[1]
	if(length(locations)<2) {
		Ny = round(locations*diff(data@bbox[2,])/diff(data@bbox[1,]))
	} else {
		Ny = locations[2]
	} 
	myExtent = 	extent(data)
	myExtent@ymax = myExtent@ymin + Ny * diff(data@bbox[1,])/Nx
	locations = raster(myExtent, Ny, Nx,
			data@proj4string)	
}
if(nrow(locations) * ncol(locations) > 10^7) warning("there are lots of cells in the prediction raster,\n this might take a very long time")

	xseq = c(xmin(locations)+xres(locations)/2, 
			xmax(locations)-xres(locations)/2, xres(locations))
	yseq = c(ymin(locations)+yres(locations)/2,
			ymax(locations)-yres(locations)/2, yres(locations))


	# the model
	requiredParams = c("variance","range","rough")
	if(!all(requiredParams %in% names(param)))
		warning("param has names", paste(names(param),collapse=","), 
				" must have ", paste(requiredParams, collapse=","))

	if(length(ycol)==1 & length(data) > 1) {
		ycol = data@data[,ycol]
	}
	
	param = fillParam(param)
	model = modelRandomFields(param)
	
	if(param["nugget"] > 0) {
		err.model = "nugget"
		err.param=c(1,param["nugget"],0)
		nuggetSd = sqrt(param["nugget"])
	} else {
		err.model = err.param = NULL
		nugget.in.prediction=FALSE
	}

	
simFun = function(D) {
	res = CondSimu(krige.method="O",
			x=xseq, y = yseq, grid=TRUE, gridtriple=TRUE,
			param=NULL, model=model,
			given=data@coords,
			data=ycol, 
			err.model=err.model,
			err.param=err.param, method="direct decomp."
	)		
	if(nugget.in.prediction){
		res= res + rnorm(length(res), sd=nuggetSd)
	}		
	values(locations) = as.vector((res[,seq(dim(res)[2], 1)]))
	if(!is.null(fun)) {
		locations = fun(locations)
	}
	locations
}		

	result = mapply(simFun, 1:Nsim)

	if(is.null(fun)) {
		resultRaster = result[[1]]
		if(Nsim > 1) {
			for(Dsim in 1:Nsim){
				resultRaster = stack(resultRaster, result[[Dsim]])
			}
		} 
		result = resultRaster			
	}
	result	
}


