stackRasterList = function(x, template=x[[1]],method='ngb') {

	if(class(x)=="SpatialPolygonsDataFrame")
		x = list(x)
	
	if(class(x)=="RasterLayer") {
		x = list(x)
		names(x) =names(x[[1]])
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
	}
 
			
	result = template
	template = as(template, "BasicRaster")
	template2 = raster(template)
	
	modefun = function(qq, na.rm=NULL) c(as.numeric(names(which.max(table(qq)))), NA)[1]
	funList = list(ngb=modefun,
			bilinear=mean)
	
	for(D in 1:Nlayers) {
 		if(class(x[[D]])=="SpatialPolygonsDataFrame"){
			if(length(names(x[[D]]))!=1)
				warning("polygon ", D, "has more than one data column, using the first" )
			result = addLayer(result,
					rasterize(
							spTransform(x[[D]], template@crs), 
							raster(template), field=names(x[[D]][1])
					)
			)
		} else { # not a spdf
		if(as(x[[D]], 'BasicRaster')==template) {
			# same projection, same resolution
			result = addLayer(result, x[[D]])			
		}	 else {
			# check to see if it's a categorical variable
			if(x[[D]]@data@isfactor) {
				method[D] = "ngb"
			} 
			
			# same projection, different resolution
			testcrs =CRS(template@crs@projargs)@projargs == CRS(x[[D]]@crs@projargs)@projargs
			if(is.na(testcrs)) testcrs = T
			if(testcrs) {
				# should we aggregate?
				toAgg = floor(min((dim(x[[D]])/dim(template2))[1:2]))
				if(toAgg > 1) {
					aggFun = funList[[method[D]]]
					thelevels = levels(x[[D]])
					x[[D]] = aggregate(x[[D]], fact=toAgg,
							fun=aggFun)
					if(!is.null(thelevels))
						levels(x[[D]]) = thelevels
				}
				
				toAdd = raster::resample(x[[D]], template2, method=method[D])
				if(!is.null(levels(x[[D]]))) {
					levels(toAdd) = levels(x[[D]])
				}
				result = addLayer(result,	toAdd)
				
			} else {
				# same resolution
			result = addLayer(result,
					projectRaster(x[[D]], template,method=method[D])
			)
		}
		}
	}
	}
	if(Nlayers == (dim(result)[3]-1) )
		result = result[[-1]]
	names(result) = names(x)
	result
}