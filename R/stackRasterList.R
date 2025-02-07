stackRasterList = function(x, template=x[[1]],method='near',mc.cores=NULL) {
# method mean and mode?  optino for a function?  defalut to 'auto'

# TO DO: fix with mc.cores > 1, currently pointer error.

	if(any(class(x)=="SpatVector"))
		x = list(x)
	
	if(any(class(x)=="SpatRaster")) {
		x = as.list(x)
		names(x) =unlist(lapply(x, names))
	}
	# TO DO: if x is list of rasters where some are multi-layered
	if(is.null(names(method))) {
	  method = rep_len(method, length(x))
	  names(method) = names(x)
	}
	if(is.list(x)) {
		if(is.null(names(x)))
			names(x) = paste("c", seq(1, length(x)),sep="")	
	}

	
	Nlayers = length(names(x))
	
	if(length(method)==Nlayers) {
		if(length(names(x)) & all( names(x)%in%names(method)))
			method = method[names(x)]
	} else {
		method = rep(method, Nlayers)
		names(method) = names(x)
	}
 	
	modefun = function(qq, na.rm=NULL) {
			# use DescTools::Mode(res)
		res = as.data.frame(table(qq))
		if(nrow(res)) {
			res = as.numeric(as.character(res[which.max(res[,2]),1]))
		} else {
			res = NA
		}
		res
	}

	funList = list(near=modefun, bilinear=mean)
	

	template = rast(template)
	template2 = rast(template)
	

	# function to reproject rasters
	projfun = function(D) {

		if(!nchar(crs(x[[D]]))) terra::crs(x[[D]]) = crs(template)
		if(any(class(x[[D]])=="SpatVector")){
			if(length(names(x[[D]]))!=1)
				warning("polygon ", D, "has more than one data column, using the first" )
			


			toAdd =  
					rasterize(
							project(x[[D]][,1], 
								crs(template)), 
							rast(template))
 			if(is.numeric(values(x[[D]])[,1])) {
# 				toAdd = deratify(toAdd)
 				toAdd = as.numeric(toAdd, 2)
 			}
		} else { # not a spatvect, assume it's a raster
			if(compareGeom(rast(x[[D]]), template, stopOnError=FALSE)) {
				# same projection, same resolution
				toAdd =  x[[D]]			
			} else { # different projection or resolution
				# check to see if it's a categorical variable
				if(any(is.factor(x[[D]]))) {
					method[D] = "near"
				} 

				thelevels = levels(x[[D]])[[1]]
				
				# same projection, different resolution
				testcrs =compareGeom(template, x[[D]],
					ext=FALSE,rowcol=FALSE,crs=TRUE,stopOnError=FALSE)				
				if(is.na(testcrs)) testcrs = TRUE
				if(testcrs) { # same resolution
					# should we aggregate?
					toAgg = floor(min(
						dim(x[[D]])[1:2]/dim(template2)[1:2]
					))
					if(toAgg > 1) {
						aggFun = funList[[method[D]]]
						if(is.factor(x[[D]])) {
							xToAgg = as.numeric(x[[D]])
						} else {
							xToAgg = x[[D]]
						}
						xagg = terra::aggregate(xToAgg, fact=toAgg,
								fun=aggFun)
						levels(xagg) = levels(x[[D]])

					} else {
						xagg = x[[D]]
					}
					toAdd = terra::resample(xagg, template2, method=method[D])
				} else { # differenet resolution
					# different resolution
					toAdd = project(x[[D]], template, method=method[D])
				} # end different resolution
				
				if(!is.null(thelevels)) {
					if(!identical(thelevels, '')) {
						levels(toAdd) = thelevels
					}
				}
			} # end different projection or resolution
		} # end not SPDF
		if(nlyr(toAdd)==1) {
			names(toAdd) = D
		} else {
			names(toAdd) = names(x[[D]])
		}
		toAdd
	} # end projfun
	
	# reproject all the rasters
	if(!is.null(mc.cores)) {
		resultList = parallel::mcmapply(
				projfun, D=names(x),
				mc.cores=mc.cores)
	} else {
		resultList = mapply(projfun, D=names(x))
	}


result = rast(resultList[1:Nlayers])	

	if(Nlayers == (dim(result)[3]-1) )
		result = result[[-1]]
	result
}