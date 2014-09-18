setClass('RMmodel', 
		representation(
				# call='RMexp(var=1, sclae=1, Aniso=id, proj=id)
				call = "language",
				# name='RMexp'
				name = "character",
				# submodels=NULL, submodels=list(RMmodel1, RMmodel2)
				submodels = "list",
				# model specific parameter 
				par.model = "list",
				# var=1, scale=1, Aniso=id, proj=id 
				par.general = "list"
		)
)


setGeneric('RFsimulate', function(model,x, data=NULL, err.model=NULL, n=1, ...) 
			standardGeneric("RFsimulate")
)


setMethod("RFsimulate", 
		signature("RMmodel", "GridTopology"),
		function(model, x,  data = NULL, 
		err.model=NULL, n = 1, ...)  {


	if (requireNamespace("RandomFields", quietly = TRUE)) { 
		
	# convert data to an RFspdf (it might be a vanilla spdf)	
	if(!is.null(data)) {
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
	
	theArgs = list(...)
	theArgs$model = model
	theArgs$x = x
	if(!is.null(data))
		theArgs$data = data
	if(!is.null(err.model))
		theArgs$err.model = err.model
	theArgs$n = n
	theArgs$spConform=TRUE
	
	res= try(do.call(RandomFields::RFsimulate, theArgs))

} else {
	warning("RandomFields package not available")
	res = NULL
}
res
}
)

setMethod("RFsimulate", 
		signature("RMmodel", "SpatialPoints"),
		function(model, x, 	data = NULL, 
				err.model=NULL, n = 1, ...)  {
			
		if (requireNamespace("RandomFields", quietly = TRUE)) { 
			# convert data to an RFspdf (it might be a vanilla spdf)	
				
			if(!is.null(data)) {
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
			}
		}

		theArgs = list(...)
		theArgs$model = model
		theArgs$x = x@coords[,1]
		theArgs$y = x@coords[,2]
		if(!is.null(data))
			theArgs$data = data
		if(!is.null(err.model))
			theArgs$err.model = err.model
		theArgs$n = n
		theArgs$spConform=TRUE
		
		res= try(do.call(RandomFields::RFsimulate, theArgs))

	} else { # RandomFields not available
		warning("RandomFields package not available")
		res = NULL
	}
	res
}
)


setMethod("RFsimulate", 
		signature("numeric", "SpatialPoints"), 
	function(model, x,  data = NULL, 
		 err.model=NULL, n = 1, ...)  {

if (requireNamespace("RandomFields", quietly = TRUE)) { 
	model = modelRandomFields(model)
	if(!is.null(err.model))
		err.model = RandomFields::RMnugget(var=err.model)
	
	theSim=callGeneric(model, x,  
		data=data,	err.model= err.model, n=n  ,  ...)
	if(class(theSim)%in%c('try-error', 'NULL')) {
		warning("error in RandomFields")
		theSim=as.data.frame(matrix(NA, length(x), n))
	} else {
	theSim = theSim@data
	}

} else { #RandomFields not available
 
		theCov = matern(x, param=model)
		if(!is.null(data)) {
			covd = matern(data, param=model)
			covpreddata = matern(x, y=data, param=model)
			if(!is.null(err.model))
				diag(covd) = diag(covd) + err.model
			Linv = solve(chol(covd))
			xcov =  tcrossprod( covpreddata,Linv)
			theCov  =theCov - tcrossprod(xcov)
		}
		theChol = chol(theCov)
		theRandom = matrix(rnorm(n*nrow(theCov)), nrow=nrow(theCov), ncol=n)
		theSim = theChol %*% theRandom
		if(!is.null(data)) {
			theSim = theSim + xcov %*% 
					(Linv %*% data.frame(data)[,1])	
		}
		theSim = as.data.frame(theSim)
	} # end no RandomFields	
	
	names(theSim) = paste("sim", 1:n, sep="")
	
	res = SpatialPointsDataFrame(SpatialPoints(x),
			data=theSim,
			proj4string=CRS(projection(x))
	)
	
	res
	}
)

setMethod("RFsimulate", signature("numeric", "GridTopology"), 
	function(model, x, data = NULL, 
				err.model=NULL, n = 1, ...)  {
		if (requireNamespace("RandomFields", quietly = TRUE)) { 
			model = modelRandomFields(model)
			if(!is.null(err.model))
				err.model = RandomFields::RMnugget(var=err.model)

			theSim=callGeneric(model, x,  
					data=data,	err.model= err.model,
					n=n  , 
					...)
			if(class(theSim)%in%c('try-error', 'NULL')) {
					warning("error in RandomFields")
					theSim=as.data.frame(matrix(NA, prod(x@cells.dim), n))
			} else {
					theSim = theSim@data
			}
		} else { #RandomFields not available
			res2 = raster(x)
			if(n>1) {
					res2 = brick(res2, nl=n)
			}
			theCov = matern(res2, param=model)
			if(!is.null(data)) {
				 covd = matern(data, param=model)
				 covpreddata = matern(res2, y=data, param=model)
				 if(!is.null(err.model))
					 diag(covd) = diag(covd) + err.model
				 Linv = solve(chol(covd))
				 xcov =  tcrossprod( covpreddata,Linv)
				 theCov  =theCov - tcrossprod(xcov)
			}
			theChol = chol(theCov)
			theRandom = matrix(rnorm(n*nrow(theCov)), nrow=nrow(theCov), 
							ncol=n)
			theSim = crossprod(theChol , theRandom)
			if(!is.null(data)) {
				theSim = theSim + xcov %*% 
							(Linv %*% data.frame(data)[,1])	
			}
			theSim = as.data.frame(theSim)

		}
		names(theSim) = paste("sim", 1:n, sep="")
		res = SpatialGridDataFrame(x,theSim)
			
		res
		}
)
			
RFsimulate.Raster = function(
		model, x,
		data=NULL, 
		err.model=NULL, n = 1, ...)  {
			
		theproj = projection(x)
		x = as(x, "GridTopology")
		res2 = callGeneric( 
				model, x,  
				data=data, err.model=err.model, n=n  , ... 
		)
		res2 = raster(res2)			
		proj4string(res2) = CRS(theproj)
			
		return(res2)
}


RFsimulate.SPgrid	=	function(
		model, x, 
				data=NULL, 
				 err.model=NULL, n = 1, ...)  {
			xOrig = x
			x= as(x, "GridTopology")

			res2 = callGeneric( 
					model, x,  
					data=data,	
					err.model= err.model, n=n  , ... 
			)
			
			if(!length(grep("DataFrame$", class(xOrig)))) {
				xOrig = as(xOrig, paste(class(xOrig), "DataFrame",sep=""))
			}
			xOrig@data = res2@data
			
			return(xOrig)
}

setMethod("RFsimulate", signature("numeric", "Raster"), 
		RFsimulate.Raster
)

setMethod("RFsimulate", 
		signature("numeric", "SpatialPixels"), 
		RFsimulate.SPgrid)


setMethod("RFsimulate", 
		signature("numeric", "SpatialGrid"), 
		RFsimulate.SPgrid)


setMethod("RFsimulate", 
		signature("data.frame", "ANY"),
		function(model, x, data=NULL, 
				err.model=NULL, n = nrow(model), ...)  {
			
		model = as(model, "matrix")
		
		callGeneric( 
				model, x, 
				data=data,	err.model= err.model,
				n=n  , ... 
		)
	}
)

setMethod("RFsimulate", 
		signature("matrix", "Spatial"),
	function(model, x,  data=NULL, 
	 err.model=NULL, n = nrow(model), ...)  {
	
 if(is.null(rownames(model)))
	 rownames(model) = paste("par", 1:nrow(model),sep="") 
 Siter = round(seq(1,nrow(model), len=n))

	if(!is.null(data)) {
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
	if(!is.null(err.model)) {
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
	result = NULL
	
	for(D in nrow(model):1) {


	resHere = 	callGeneric(
			model[D,], 
			x, data= data[,D], 
			err.model= err.model[D], n=1,
			...)
	result = cbind( 
		 resHere@data[,1], result
		)
	}
	
	
	resHere@data	= as.data.frame( result)

	if(!is.null(rownames(model)))
		names(resHere) = rownames(model)
	resHere
}
)

setMethod("RFsimulate", 
		signature("matrix", "Raster"),
		function(model, x, data=NULL, 
				 err.model=NULL, n = nrow(model), ...)  {

			 if(is.null(rownames(model)))
				rownames(model) = paste("par", 1:nrow(model),sep="") 
			Siter = round(seq(1,nrow(model), len=n))
			
			if(!is.null(data)) {
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
			if(!is.null(err.model)) {
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
			result = NULL
			
			for(D in nrow(model):1) {
				resultHere = 	callGeneric(
						model[D,], 
						x, data= data[,D], 
						err.model[D], n=1,
						...)
				result = brick( 
						resultHere,result
				)
			}
			
			if(!is.null(rownames(model)))
				names(result) = rownames(model)
			result
		}
)

