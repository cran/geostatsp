GaussRF = function(x,param=c(variance=1, range=1, rough=1), ...) {
	UseMethod("GaussRF")
	
}

GaussRF.Raster = function(x,param=c(variance=1, range=1, rough=1), ...){

	xseq = seq(x@extent@xmin, x@extent@xmax, len=x@ncols)
	yseq = seq(x@extent@ymin, x@extent@ymax, len=x@nrows)

 	
	
	
	res = GaussRF(x=xseq, param=param,
			y=yseq,	grid=TRUE, 
		...
			)

			
	resRast = raster(t(res[,seq(dim(res)[2], 1)]),
			x@extent@xmin, x@extent@xmax,
			x@extent@ymin, x@extent@ymax, crs=x@crs)
	
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

		param["scale"] = 2*param["range"] 
		
		if(any(names(param)=="aniso.ratio")){
			# geometric anisotropy
			if(any(names(param)=="aniso.angle.degrees") & 
					!any(names(param)=="aniso.angle.radians") ) {
				param["aniso.angle.radians"] = param["aniso.angle.degrees"]*2*pi/360				
		}
			model=list("$", var=param["variance"],   
					A=anisoMatrix(param["aniso.angle.radians"],
							param["scale"]*c(param["aniso.ratio"],1)
				),
				list("matern", nu=param["rough"]))	
		} else {
		
			model=list("$", var=param["variance"],   
					s=param["scale"],
					list("matern", nu=param["rough"]))	
			
		}	
		theArgs$model = model	
	} else {
		if(!is.list(theArgs$model)) # if model is a list, it's the extended model definition which doesnt need the param argument
			theArgs$param = param
	}
	
	result = do.call( RandomFields::GaussRF, theArgs)
	if(any(names(param)=="mean") )  {
		result = result + param["mean"]
	}
	result	
	
}


