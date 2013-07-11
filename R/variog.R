variog = function(geodata, ...) {
	UseMethod("variog")
	
}


variog.default <- function(geodata, ...) {
	geoR::variog(geodata, ...)
}
			
variog.SpatialPointsDataFrame = function(geodata,formula, ...) {
	
	theResid = lm(formula, data=geodata@data)$resid
 
	result = geoR::variog(coords=geodata@coords, data=theResid, ...)
 
	result
}


variog.mc.env = function(geodata, formula, ...) {
	UseMethod("variog.mc.env")
}

variog.mc.env.default = function(geodata, ...) {
	geoR::variog.mc.env(geodata,...) 
}

variog.mc.env.SpatialPointsDataFrame = function(geodata,formula, ...) {
	
	theResid = lm(formula, data=geodata@data)$resid
	
	geoR::variog.mc.env(coords=geodata@coords, data=theResid, ...)
	
}
