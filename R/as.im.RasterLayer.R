as.im = function(x, ...) {
	UseMethod("as.im")	
}

as.im.default = function(x, ...){
	spatstat::as.im(x, ...)			
}

as.im.RasterLayer = function(x, ...) {
	spatstat::as.im(as.matrix(x)[nrow(x):1,], 
			xrange=bbox(x)[1,], 
			yrange=bbox(x)[2,], ...)
}