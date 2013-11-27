
matern = function(x, y=NULL, param=c(range=1, variance=1, shape=1)) {
	UseMethod("matern")
	
}

matern.dist = function( x, y=NULL, param=c(range=1, variance=1, shape=1)) {

	param=fillParam(param)
	resultVec = matern(as.vector(x), param=param)
	resultMat = matrix(0, attributes(x)$Size, 
			attributes(x)$Size)
	resultMat[lower.tri(resultMat)] = resultVec

	result = new("dsyMatrix", 
			Dim = dim(resultMat), uplo="L",
			x=c(resultMat))
	Matrix::diag(result) = attributes(resultVec)$param["variance"]
	
	if(attributes(resultVec)$zeros/length(resultVec)>0.5) {
		result@x[result@x < result[1,1]*1e-06] = 0
		result = as(result, "dsCMatrix")
	}
	
	attributes(result)$param = attributes(resultVec)$param	
	
	result
}

matern.SpatialPointsDataFrame = function( x, y=NULL, param=c(range=1, variance=1, shape=1)) {
	x = SpatialPoints(x)
	matern(x=x, y=y, param=param)
}


matern.Raster = function( x, y=NULL, param=c(range=1, variance=1, shape=1))
 {
	param = fillParam(param)
	 if(is.null(y)) {
		 y=x
		 symm=TRUE
	 } else {
		 symm=FALSE
	 }
	 # convert  y to spatial points, no matter what it is
	 if(is.vector(y)) y = matrix(y[1:2], 1,2) 
	 y = SpatialPoints(y)

	 Ny = length(y)
	 
	 
	 resC= .C("maternArasterBpoints", 
			 as.double(xmin(x)), as.double(xres(x)), as.integer(ncol(x)), 
			 as.double(ymax(x)), as.double(yres(x)), as.integer(nrow(x)),
			 as.double(y@coords[,1]), as.double(y@coords[,2]), 
			 N=as.integer(Ny), 
			 result=as.double(array(0, c(nrow(x),ncol(x),Ny))),
			 xscale=as.double(param["range"]),
			 varscale=as.double(param["shape"]),
			 as.double(param["variance"]),
			 as.double(param["anisoRatio"]),
			 as.double(param["anisoAngleRadians"])
	 )
	
	if(Ny ==1) {
		values(x) = resC$result		
	} else {
		x = matrix(resC$result, nrow=ncell(x), ncol=Ny)
		if(symm){
			if(resC$N / length(x) > 0.5) {
				x[x < param["variance"]*1e-06] = 0
				x = as(x, "dsCMatrix")
			} else {
			x = as(x, "dsyMatrix")
			}
		} else {
			if(resC$N / length(x) > 0.5) {
				# convert to sparse matrix
				x[x < param["variance"]*1e-06] = 0
				x = as(x, "dgCMatrix")
			}
		}
	} 
	attributes(x)$param = param	 
	x

}



matern.SpatialPoints = function(x, y=NULL,param=c(range=1, variance=1, shape=1)
		){

	param = fillParam(param)		
			
	if(!is.null(y)) {	
		# haven't written this in C yet.. rotate and create distances in R
		if(length(grep("SpatialPoints", class(y)))) {
			y = y@coords[,1] + 1i*y@coords[,2]  
		}
		if(length(grep("^Raster", class(y)))) {
			y = as.data.frame(y, xy=TRUE)
			y = y[,"x"] + 1i*y[,"y"]  
		}
		
		if(length(y)==2 & !is.complex(y)){
			y = y[1] + 1i*y[2]
		}

		x = x@coords[,1] + 1i*x@coords[,2]
		
			
		x = x * exp(1i*param["anisoAngleRadians"])
		x = Re(x) +  (1i/ param["anisoRatio"] )*Im(x)
		y = y * exp(1i*param["anisoAngleRadians"])
		y = Re(y) +  (1i/ param["anisoRatio"] )*Im(y)
				
		thedist = Mod(outer(x, y, FUN="-"))
		result= matern(x=thedist, y=NULL, param=param)

 		
	} else { # y is null
#	void maternAniso(double *x, double *y, long *N,
#					double *result,
#					double  *range, double*shape, 
#	double *variance,
#				double *anisoRatio, double *anisoAngleRadians) {
					
 	resC = .C("maternAniso", 
			as.double(x@coords[,1]),
				as.double(x@coords[,2]), 
				N= as.integer(length(x)),
				result=as.double(rep(-99.9, length(x)^2)),
				as.double(param["range"]),
				as.double(param["shape"]),
				as.double(param["variance"]),
				as.double(param["anisoRatio"]),
				as.double(param["anisoAngleRadians"])
			)
	result = new("dsyMatrix", 
				Dim = c(length(x), length(x)), uplo="L",
						x=resC$result)
	Matrix::diag(result) = param["variance"]
				
	if(resC$N/length(x)^2>0.25) {
			result@x[result@x < result[1,1]*1e-06] = 0
			result = as(result, "dsCMatrix")
	}
				
				
	attributes(result)$param = param		
	}
		
	result
}

matern.default = function( x, y=NULL,param=c(range=1, variance=1, shape=1))
{
	# x is distances (matrix or vector), y is ignored	
	names(param) = gsub("^var$", "variance", names(param))
	
	if(!any(names(param)=="variance") & any(names(param)=="sdSpatial"))
		param["variance"]= param["sdSpatial"]^2
	
	haveVariance = any(names(param)=="variance")
	if(!haveVariance) 
		param["variance"]=1

	if(is.data.frame(x))
		x = as.matrix(x)	
#	void matern(double *distance, long *N,
#					double *range, double *shape, double *variance) {
	resultFull = .C("matern", as.double(x), as.integer(length(x)),
			as.double(param["range"]), as.double(param["shape"]),
			as.double(param["variance"]))
	result = resultFull[[1]]
	if(is.matrix(x)) 
		result = matrix(result, ncol=ncol(x), nrow=nrow(x))
	
	attributes(result)$param = param
	attributes(result)$zeros = resultFull[[2]]
	
	result
	
}



oldmatern = function( x, param=c(range=1, variance=1, shape=1))
{
	# R code instead of C
	# x is distances (matrix or vector), y is ignored, assume isotropic	
	param = fillParam(param)
	
	if(is.data.frame(x))
		x = as.matrix(x)	
	# do this bit in C?
	xscale = abs(x)*(sqrt(8*param["shape"])/ param["range"])
	result = ( param["variance"]/(gamma(param["shape"])* 2^(param["shape"]-1)  ) ) * 
			( xscale^param["shape"] *
				besselK(xscale , param["shape"]) )
	result[xscale==0] = 
			param["variance"] 
	result[xscale==Inf] = 0
	
	
	
	attributes(result)$param = param
	attributes(result)$xscale = (sqrt(8*param["shape"])/ param["range"])
	attributes(result)$varscale =  param["variance"]/(gamma(param["shape"])* 
					2^(param["shape"]-1)  )
	result
	
}

