matern = function(x, y=NULL, param=c(range=1, variance=1, rough=1)) {
	UseMethod("matern")
	
}

matern.dist = function( x, y=NULL, param=c(range=1, variance=1, rough=1)) {

	resultVec = matern(as.vector(x), param=param)
	resultMat = matrix(0, attributes(x)$Size, 
			attributes(x)$Size)
	resultMat[lower.tri(resultMat)] = resultVec

	result = new("dsyMatrix", 
			Dim = dim(resultMat), uplo="L",
			x=c(resultMat))
	Matrix::diag(result) = attributes(resultVec)$param["variance"]
	attributes(result)$param = attributes(resultVec)$param	
	result
}

matern.SpatialPointsDataFrame = function( x, y=NULL, param=c(range=1, variance=1, rough=1)) {
	x = SpatialPoints(x)
	matern(x=x, y=y, param=param)
}


matern.Raster = function( x, y=NULL, param=c(range=1, variance=1, rough=1))
 {

	 if(length(grep("^Raster", class(y)))){
		 y = SpatialPoints(as.data.frame(y, xy=TRUE)[,c("x","y")])
	 }
	 
	resRast = x
	x= SpatialPoints(as.data.frame(x, xy=TRUE)[,c("x","y")])

	# for some reason
	# NextMethod("matern")
	# doesn't work, it calls matern.default
	result= matern(x=x, y=y, param=param)

	if(length(result)==ncell(resRast)) {
			values(resRast) = result 
			result=resRast
	} 
	result
}



matern.SpatialPoints = function(x, y=NULL,param=c(range=1, variance=1, rough=1)
		){
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
		
		if(any(names(param)=="aniso.ratio") ) { 
			
			# geometric anisotropy
		# rotate coordinates
			if(any(names(param)=="aniso.angle.degrees") & 
				!any(names(param)=="aniso.angle.radians") ) {
			param["aniso.angle.radians"] = param["aniso.angle.degrees"]*2*pi/360				
			}
			if(!any(names(param)=="aniso.angle.radians") )
				warning("anisotropy angle not supplied")
			
			x = x * exp(1i*param["aniso.angle.radians"])
			x = Re(x) +  (1i/ param["aniso.ratio"] )*Im(x)
			if(!is.null(y)) {
				y = y * exp(1i*param["aniso.angle.radians"])
				y = Re(y) +  (1i/ param["aniso.ratio"] )*Im(y)
			}
 
		} 

  		if(!is.null(y)) {
			thedist = Mod(outer(x, y, FUN="-"))
			result= matern(x=thedist, y=NULL, param=param)

		} else {

	
			thedist = dist(cbind(Re(x), Im(x)))
		
			result = matern(x=thedist, y=NULL, param=param)
 		
		}
				
		result
}

matern.default = function( x, y=NULL,param=c(range=1, variance=1, rough=1))
{
	# x is distances, y is ignored	
	names(param) = gsub("^var$", "variance", names(param))
	
	if(!any(names(param)=="variance") & any(names(param)=="sdSpatial"))
		param["variance"]= param["sdSpatial"]^2
	

	
	haveVariance = any(names(param)=="variance")

	if(!haveVariance) 
		param["variance"]=1
	
		
	xscale = abs(x)*(sqrt(8*param["rough"])/ param["range"])
	result = ( param["variance"]/(gamma(param["rough"])* 2^(param["rough"]-1)  ) ) * 
			( xscale^param["rough"] *
				besselK(xscale , param["rough"]) )

	result[xscale==0] = 
			param["variance"] 
	
	result[xscale==Inf] = 0
	
	attributes(result)$param = param
	result
	
}
