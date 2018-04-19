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
	
	callGeneric(x = extent(x), cells, buffer)
})

setMethod("squareRaster", 
  signature(x="Extent"),
  function(x, cells=NULL, buffer=0) {
	
	if(is.null(cells))
		warning("cells must be specified if x is a matrix or extent")
	cells = as.integer(cells[1])
	x = raster(x, ncol=cells, nrow=cells)
	callGeneric(x, cells, buffer)
})

setMethod("squareRaster", 
  signature(x="Raster"), 
  function(x, cells=NULL, buffer=0) {
	x=raster(x)
	if(is.null(cells)) {
		cells = ncol(x)
	} else {
		ncol(x) = as.integer(cells[1])
	}
	Ny = ceiling(signif( (ymax(x) - ymin(x))/xres(x), 10) )
	ymax(x) = ymin(x) + Ny*xres(x)
	nrow(x) = Ny
	extend(x, round(buffer/xres(x)))
})

setMethod("squareRaster", 
  signature(x="Spatial"),
		function(x, cells=NULL, buffer=0) {
	
	result = squareRaster(extent(x), cells, buffer)
	proj4string(result) = proj4string(x)
	result
})