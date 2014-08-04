
RFsimulate = function(model, x, y = NULL, z = NULL, T = NULL, grid, data, 
		distances, dim, err.model, n = 1, ...) {
	UseMethod("RFsimulate")
}


RFsimulate.RMmodel = 
		function (model, x, y = NULL, z = NULL, T = NULL, grid, data, 
				distances, dim, err.model, n = 1, ...) 
{

	# if x has a projection, copy it to the output
	theProj = try(proj4string(x), silent=TRUE)
	if(class(theProj)=="try-error")
		theProj=NA
	
	xOrig=x
	if(any(attributes(class(x))$package=="raster")) {
		x = as(x, "GridTopology")
		isRaster = TRUE
	} else {
		isRaster=FALSE
	}
	
	if(length(grep("^SpatialPoints", class(x)))){
		y=coordinates(x)[,2]
		x=coordinates(x)[,1]
		isSpatialPoints = TRUE
	} else {
		isSpatialPoints = FALSE
	}

	if(length(grep("^Spatial(Pixels|Grid)", class(x)))){
		x = getGridTopology(x)
	}
	
	# convert data to an RFspdf (it might be a vanilla spdf)	
	if(!missing(data)) {
	if(class(data)=='SpatialPointsDataFrame'){
		data = RandomFields::conventional2RFspDataFrame(
				array(
						as.matrix(data@data),
						c(dim(data@data),1)
				),
				coordinates(data),
				n=dim(data)[2]
				)
#		data=as(data, "RFspatialPointsDataFrame") 
	} }


	res2=RandomFields::RFsimulate(
			model, x, y, z , T , grid, data, 
			distances, dim, err.model, n , ...
			)
		
	# assign a proj4string if necessary
	if(class(try(proj4string(res2),silent=TRUE))!="try-error") {
		if(is.na(proj4string(res2))) {
			if(class(try(CRS(theProj), silent=TRUE))!= "try-error") {
					proj4string(res2) = CRS(theProj)
			}
		}				
	}

	
	# if x is a raster and the raster package is available
	# convert the results to a raster
	if(isRaster & 
			any(rownames(installed.packages())=="raster")) {

		if(n==1) {
			# one simulation, convert to raster
			if(is.matrix(res2)) {
				res3 = xOrig
				values(res3) = as.vector(res2)
			} else  {
				res3 = raster::raster(res2)
			}
		} else {
			# more than one simulate, convert to raster brick
			res3 = list()
			for(D in 1:n) {
				res3[[D]] = 
						raster::raster(res2, layer=D)
			}
			res3= do.call(raster::brick, res3)
		}
		res2 = res3
	} 
	res2
}

RFsimulate.numeric = function(model, x, y = NULL, z = NULL, T = NULL, grid, data, 
		distances, dim, err.model, n = 1, ...)  {

	model = modelRandomFields(model)

	# if err.model is numeric, assume it's a nugget variance
	if(!missing(err.model)){
		if(is.numeric(err.model)) {
			err.model = 
					RandomFields::RMnugget(
							var=err.model)
		}
	}
	
	
	if(!missing(data)) {
		if(class(data)=="SpatialPointsDataFrame") {
			if(ncol(data)==1){
				data = RandomFields::conventional2RFspDataFrame(
						data.frame(data)[,1],
						coords=coordinates(data),
						n=1, vdim=1
				) 
			}
		}
	}
	
	
	#RandomFields::
	#
	res2=geostatsp::RFsimulate(model, x, y  , z  , T   , grid, 
		data,	distances, dim, err.model, n  , ...)

# convert to non-RandomFields object
if(length(grep("[Ss]patialPointsDataFrame", class(res2)))){
	res2 = as(res2, "SpatialPointsDataFrame")
}

	return(res2)

	
}

RFsimulate.data.frame = function(model, x, y = NULL, z = NULL, T = NULL, grid, data, 
		distances, dim, err.model, n = 1, ...)  {
	model = as(model, "matrix")
	geostatsp::RFsimulate(model, x, y  , z  , T   , grid, data, 
			distances, dim, err.model, n  , ...)
}

RFsimulate.matrix = function(model, x, y = NULL, z = NULL, T = NULL, grid, data, 
		distances, dim, err.model, n = 1, ...)  {
	

	Siter = round(seq(1,nrow(model), len=n))

	if(!missing(data)) {
		if(class(data)!= "SpatialPointsDataFrame")
			warning("data should be a SpatialPointsDataFrame")
	# check data variables
	if(ncol(data) == 1) {
		data = data[,rep(1,length(Siter))]
		
	} else if(ncol(data) > 1) {
	# if there's more than one, assume we're interested in the first one
	
		if(ncol(data)!= nrow(model)){
			warning("number of columnns in data should be either 1 or equal to number of rows of model")
		}
		data = data[,Siter]		
	}
	} else {
		# do something so data[,D] doesn't break
		data = NULL
	}
	if(!missing(err.model)) {
		if(is.numeric(err.model)){
			if(length(err.model)==1) {
				err.model = rep(err.model, length(Siter))
			} else if(length(err.model)==nrow(model)){
				err.model = err.model[Siter]
			} else {
				warning("number of values in err.model should be either 1 or equal to number of rows of model")
			}
		}
	} else {
		err.model= NULL
	}

	model= model[Siter,]
	
	result = list()
	
	
	for(D in 1:nrow(model)) {

	result[[D]] = geostatsp::RFsimulate(
			model[D,], 
			x, y  , z  , T   , grid, data[,D], 
			distances, dim, err.model[D], n=1  , ...)

	}
	
	# try to convert result from a list to raster brick or spdf
	if(attributes(class(result[[1]]))$package=="raster") {
		result = do.call(raster::brick, result)
	}
	if(length(grep("Spatial*DataFrame", class(result[[1]]) )))	{
		result[[1]]@data	= as.data.frame(
				lapply(result, function(qq) data.frame(qq)[,1])
		)
		result = result[[1]]
	}
	
	result
}

