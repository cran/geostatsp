setGeneric('squareRaster', 
  function(
    x, cells=NULL, buffer=0) {
    standardGeneric("squareRaster")
  }
  )



setMethod("squareRaster", 
  signature(x="matrix"), 
  function(x,   cells=NULL, buffer=0) {
# assume it's a bbox	
	
	callGeneric(x = ext(x), cells, buffer)
})

setMethod("squareRaster", 
  signature(x="SpatExtent"),
  function(x, cells=NULL, buffer=0) {
	
	if(is.null(cells))
		warning("cells must be specified if x is a matrix or extent")
	cells = as.integer(cells[1])
	x = rast(extent=x, ncol=cells, nrow=cells)
	callGeneric(x, cells, buffer)
})

setMethod("squareRaster", 
  signature(x="SpatRaster"), 
  function(x, cells=NULL, buffer=0) {
	x=rast(x, nl=1)
	if(is.null(cells)) {
		cells = ncol(x)
	} else {
		terra::ncol(x) = as.integer(cells[1])
	}
	Ny = ceiling(signif( (ymax(x) - ymin(x))/xres(x), 10) )
	terra::ymax(x) = ymin(x) + Ny*xres(x)
	terra::nrow(x) = Ny
	extend(x, ceiling(buffer/xres(x)))
})

setMethod("squareRaster", 
  signature(x="SpatVector"),
		function(x, cells=NULL, buffer=0) {
	
	result = squareRaster(ext(x), cells, buffer)
#	proj4string(result) = proj4string(x)
	terra::crs(result) = crs(x)
	result
})