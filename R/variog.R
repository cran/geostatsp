variog = function(geodata, ...) {
	UseMethod("variog")	
}

variog.default <- function(geodata, ...) {
	if (requireNamespace("geoR", quietly = TRUE)) { 
		result=geoR::variog(geodata, ...) 
	} else {
		result = NULL
	}
result
}
			
variog.SpatialPointsDataFrame = function(geodata,formula, ...) {
	
	theResid = lm(formula, data=geodata@data)$resid
 
	if (requireNamespace("geoR", quietly = TRUE)) { 
		
		result = geoR::variog(coords=geodata@coords, data=theResid, ...)
 
	} else {
		result = NULL
	}
	result
}


variog.mc.env = function(geodata, formula, ...) {
	UseMethod("variog.mc.env")
}

variog.mc.env.default = function(geodata, ...) {
	if (requireNamespace("geoR", quietly = TRUE)) { 
		result = geoR::variog.mc.env(geodata,...) 
	} else {
		result = NULL
	}
	result
}

variog.mc.env.SpatialPointsDataFrame = function(geodata,formula, ...) {
	
	theResid = lm(formula, data=geodata@data)$resid
	if (requireNamespace("geoR", quietly = TRUE)) { 
		result = geoR::variog.mc.env(coords=geodata@coords, data=theResid, ...)
	} else {
			result = NULL
	}
result
}