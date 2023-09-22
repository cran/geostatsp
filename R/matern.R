typeIntFromString = function(type) {
	if(is.character(type)) {
		type = gsub("iance$|esky$|ision", "", tolower(type)[1])    
		type = as.integer(c(var=1,chol=2,prec=3,inversechol=4)[type])    
	}
	type
}


matern = function(x, 
	param=c(range=1, variance=1, shape=1), 
	type=c('variance','cholesky','precision','inverseCholesky'),y=NULL) {
	UseMethod("matern")
}

matern.dist = function(x,
	param=c(range=1, variance=1, shape=1),
	type=c('variance','cholesky','precision','inverseCholesky'), y=NULL) {
	
	type = typeIntFromString(type)


	param=fillParam(param)
	x = as.matrix(x)
	cres = .C(C_matern, 
		x=as.double(x), 
		N=as.integer(ncol(x)),
		result = as.double(rep(-99.9, prod(dim(x)))),
		as.double(param["range"]), 
		shape=as.double(param["shape"]),
		as.double(param["variance"]),
		as.double(param["nugget"]),
		type=as.integer(type),
		halfLogDet=as.double(-9.9))

	
	if (type==2 | type==4){
		result = new("dtrMatrix", 
				Dim = dim(x), 
				uplo="L",
				x=cres$result)
		
		attributes(result)$logDetHalf = cres$halfLogDet	
		attributes(result)$cholInfo = cres$type

	} else {
		result = new("dpoMatrix", 
			Dim = dim(x), 
			uplo="L",
			x=cres$result)
	}

	attributes(result)$param = param	
	result
}

matern.dsyMatrix = function(x, 
	param=c(range=1, variance=1, shape=1),
	type=c('variance','cholesky','precision','inverseCholesky'),
	y=NULL) {
	
	param=fillParam(param)[c('range','shape','variance','nugget')]

	type = typeIntFromString(type)


	N = nrow(x)
	if(type %in% c(2,4)) { # triangular matrix
		result = new('dtrMatrix', uplo = 'L', diag = 'N', x = rep(0.0, N*N), Dim = rep(N,2))
	} else {
		result = new('dsyMatrix', uplo = 'L', x = rep(0.0, N*N), Dim = rep(N,2))
	}
	dimnames(result) = dimnames(x)

	halfLogDet = .Call(
		C_maternDistance,
		x, result, param, type
		)

	attributes(result)$param = param
	attributes(result)$type = type
	attributes(result)$halfLogDet = halfLogDet

	result	
}


matern.SpatVector = function(x, 
	param=c(range=1, variance=1, shape=1), 
	type=c('variance','cholesky','precision','inverseCholesky'), 
	y=NULL) {

	typeOrig = type
	type = typeIntFromString(type)

	N = length(x)
	if(type %in% c(2,4)) { # triangular matrix
		result = new('dtrMatrix', uplo = 'L', diag = 'N', x = rep(0.0, N*N), Dim = rep(N,2))
	} else {
		result = new('dsyMatrix', uplo = 'L', x = rep(0.0, N*N), Dim = rep(N,2))
	}

	halfLogDet = .Call(C_maternPoints,
		crds(x),
		result,
		fillParam(param)[c(
			'range','shape','variance',
			'anisoRatio','anisoAngleRadians','nugget')
		],
		as.integer(type))
	if(type >1) attributes(result)$halfLogDet = halfLogDet
	attributes(result)$type = typeOrig
	attributes(result)$param = param
	result
}


matern.SpatRaster = function(x, 
	param=c(range=1, variance=1, shape=1),
	type=c('variance','cholesky','precision','inverseCholesky'), 
	y=NULL) {

	type = typeIntFromString(type)
	
	param = fillParam(param)
	if(is.null(y)) {
		suppressWarnings(y <- as.points(x))
		x = matern(y, param, type)
	} else {
		# convert  y to spatial points, no matter what it is
		if(is.vector(y)) y = matrix(y[1:2], 1,2)
		if(is.matrix(y)) {
			y = vect(y, crs=crs(x))
		}
		suppressWarnings(y <- as.points(y))		
	
		Ny = length(y)
	
		resC= .C(C_maternArasterBpoints, 
			as.double(xmin(x)), 
			as.double(xres(x)), 
			as.integer(ncol(x)), 
			as.double(ymax(x)),
			as.double(yres(x)), 
			as.integer(nrow(x)),
			as.double(crds(y)[,1]), 
			as.double(crds(y)[,2]), 
			N=as.integer(Ny), 
			result=as.double(array(0, c(nrow(x),ncol(x),Ny))),
			xscale=as.double(param["range"]),
			varscale=as.double(param["shape"]),
			as.double(param["variance"]),
			as.double(param["anisoRatio"]),
			as.double(param["anisoAngleRadians"])
		)
	
		if(Ny ==1) {
			terra::values(x) = resC$result		
		} else {
			x = matrix(resC$result, nrow=ncell(x), ncol=Ny)
		}
	}
	attributes(x)$param = param
	x
}



notExportedMaternSpatialPointsXX = function(x,
	param=c(range=1, variance=1, shape=1),
	type=c('variance','cholesky','precision','inverseCholesky'), y=NULL
	){
	
	type = typeIntFromString(type)

	param = fillParam(param)		
	
	if(!is.null(y)) {	
		# haven't written this in C yet.. rotate and create distances in R
		if(length(grep("SpatVector", class(y)))) {
			y = crds(y)[,1] + 1i*crds(y)[,2]  
		}
		if(length(grep("^SpatRaster", class(y)))) {
			y = crds(y)
			y = y[,"x"] + 1i*y[,"y"]  
		}
		
		if(length(y)==2 & !is.complex(y)){
			y = y[1] + 1i*y[2]
		}
		
		x = crds(x)[,1] + 1i*crds(x)[,2]
		
		
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

	resC = .C(C_maternAniso, 
		as.double(crds(x)[,1]),
		as.double(crds(x)[,2]), 
		N= as.integer(length(x)),
		result=as.double(rep(-99.9, length(x)^2)),
		as.double(param["range"]),
		shape=as.double(param["shape"]),
		as.double(param["variance"]),
		as.double(param["anisoRatio"]),
		as.double(param["anisoAngleRadians"]),
		as.double(param["nugget"]),
		type=as.integer(type),
		halfLogDet=as.double(-9.9)
		)
	if(type==2 | type==4){
		result = as(
			new("dtrMatrix", 
				Dim = as.integer(c(length(x), length(x))), 
				uplo="L",
				x=resC$result),
			"Cholesky")
		attributes(result)$logDetHalf = resC$halfLogDet
		attributes(result)$cholInfo = resC$type
	} else {
		result = new("dpoMatrix", 
			Dim = as.integer(c(length(x), length(x))), 
			uplo="L",
			x=resC$result)
		if(type==3) 		
			attributes(result)$logDetHalf = resC$halfLogDet
	}

} # end y null

attributes(result)$param = param		

result
}

matern.default = function(x, 
	param=c(range=1, variance=1, shape=1),
	type=c('variance','cholesky','precision','inverseCholesky'), y=NULL) {
	# x is distances (matrix or vector), y is ignored	
	
	if(!any(names(param)=="variance") & any(names(param)=="sdSpatial"))
		param["variance"]= param["sdSpatial"]^2

	param=fillParam(param)

	if(is.data.frame(x))
		x = as.matrix(x)	
#	void matern(double *distance, long *N,
#					double *range, double *shape, double *variance) {
	resultFull = .C(
		C_matern, 
		as.double(x), 
		as.integer(length(x)),
		result= as.double(rep(-99.9, length(x))),
		as.double(param["range"]), 
		as.double(param["shape"]),
		as.double(param["variance"]),
		as.double(0), # nugget
		type=as.integer(0),
		halfLogDet=as.double(-9.9)
		)
	result = resultFull$result
	if(is.matrix(x)) 
		result = matrix(result, ncol=ncol(x), nrow=nrow(x))
	
	attributes(result)$param = param
	
	result
}
