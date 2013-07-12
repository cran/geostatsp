GaussRF = function(x, ...) {
	UseMethod("GaussRF")
	
}

GaussRF.Raster = function(x, ...){

	xseq = seq(x@extent@xmin, x@extent@xmax, len=x@ncols)
	yseq = seq(x@extent@ymin, x@extent@ymax, len=x@nrows)

	thediffx = diff(range(diff(xseq)))
	if(thediffx > 0) {
		thediffx = signif(diff(xseq[1:2]), 10^-ceiling(log10(thediffx)+1))
		xseq = seq(x@extent@xmin, by=thediffx, len=x@ncols )
	}

	thediffy = diff(range(diff(yseq)))
	if(thediffy > 0) {
		thediffy = signif(diff(yseq[1:2]), 10^-ceiling(log10(thediffy)+1))
		yseq = seq(x@extent@ymin, by=thediffy, len=x@ncols )
	}
	
	
	
	res =RandomFields::GaussRF(
			x=xseq, y=yseq,
			grid=TRUE, 
		...
			)

			
	resRast = raster(t(res[,seq(dim(res)[2], 1)]),
			x@extent@xmin, x@extent@xmax,
			x@extent@ymin, x@extent@ymax, crs=x@crs)
	
	return(resRast)
}

GaussRF.SpatialPointsDataFrame = function(x, ...){
	
	result =RandomFields::GaussRF(
		x@coords,
		...
	)

	return(result)
}

GaussRF.SpatialPoints= function(x, ...){
	
	result =RandomFields::GaussRF(
			x@coords,
			...
	)
	
	return(result)
}

GaussRF.default = function(x,  ...){
	

	RandomFields::GaussRF(x,
			...
	)
	
	
}


