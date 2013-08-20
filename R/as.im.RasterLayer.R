
as.im.RasterLayer = function(X, ...) {
	spatstat::as.im.matrix(raster::as.matrix(X)[nrow(X):1,], 
			xrange=bbox(X)[1,], 
			yrange=bbox(X)[2,], ...)
}