

GaussRF = function(x,param=c(variance=1, range=1, rough=1), ...) {
	UseMethod("GaussRF")
	
}

GaussRF.Raster = function(x,param=c(variance=1, range=1, rough=1), ...){

	xseq = c(xmin(x)+xres(x)/2, xmax(x)-xres(x)/2, xres(x))
	yseq = c(ymin(x)+yres(x)/2, ymax(x)-yres(x)/2, yres(x))
	
 	
	
	
	res = GaussRF(x=xseq, param=param,
			y=yseq,	grid=TRUE, gridtriple=TRUE,
			...
			)

			
	resRast = raster(t(res[,seq(dim(res)[2], 1)]),
			x@extent@xmin, x@extent@xmax,
			x@extent@ymin, x@extent@ymax, crs=x@crs)
	
	#attributes(resRast)$param = attributes(res)$param
	
	return(resRast)
}

GaussRF.SpatialPointsDataFrame = function(x,param=c(variance=1, range=1, rough=1), ...){
	
	x=coordinates(x)
	NextMethod("GaussRF")
}

GaussRF.SpatialPoints= function(x,param=c(variance=1, range=1, rough=1), ...){
	
	x=coordinates(x)
	NextMethod("GaussRF")
}

GaussRF.default = function(x,param=c(variance=1, range=1, rough=1),  ...){
	
	
 	theArgs = list(...)
	theArgs$x = x
	
	if(!any(names(theArgs)=="model")){
		# param is for geostatsp, not RandomFields 

		requiredParams = c("variance","range","rough")
		if(!all(requiredParams %in% names(param)))
			warning("param has names", paste(names(param),collapse=","), 
					" must have ", paste(requiredParams, collapse=","))

		
		model  = modelRandomFields(param)
		
		theArgs$model = model	
	} else {
		if(!is.list(theArgs$model)) # if model is a list, it's the extended model definition which doesnt need the param argument
			theArgs$param = param
	}
	
 
	
	result = do.call( RandomFields::GaussRF, theArgs)

	# some things break (such as spplot) if I add this as an attribute
	#attributes(result)$param = param
	
	result	
	
}


