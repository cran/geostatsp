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


variogMcEnv = function(geodata, formula, ...) {
	UseMethod("variogMcEnv")
}

variogMcEnv.default = function(geodata, ...) {
	if (requireNamespace("geoR", quietly = TRUE)) { 
		result = geoR::variog.mc.env(geodata,...) 
	} else {
		result = NULL
	}
	result
}

variogMcEnv.SpatialPointsDataFrame = function(geodata,formula, ...) {
	
	theResid = lm(formula, data=geodata@data)$resid
	if (requireNamespace("geoR", quietly = TRUE)) { 
		result = geoR::variog.mc.env(coords=geodata@coords, data=theResid, ...)
	} else {
			result = NULL
	}
result
}