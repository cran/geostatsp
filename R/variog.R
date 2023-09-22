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
			
variog.SpatVector = function(geodata,formula, ...) {
	
  rownames(terra::values(geodata)) =  1:length(geodata)
  theCoords = crds(geodata)
  rownames(theCoords) = 1:length(geodata)
  
	theResid = lm(formula, data=values(geodata))$resid
 
	if (requireNamespace("geoR", quietly = TRUE)) { 
		
		result = geoR::variog(coords=theCoords[names(theResid),], 
        data=theResid, ...)
 
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

variogMcEnv.SpatVector = function(geodata,formula, ...) {
  
  rownames(terra::values(geodata)) =  1:length(geodata)
  theCoords = crds(geodata)
  rownames(theCoords) = 1:length(geodata)
  
  
	theResid = lm(formula, data=values(geodata))$resid
	if (requireNamespace("geoR", quietly = TRUE)) { 
		result = geoR::variog.mc.env(coords=theCoords[names(theResid),], data=theResid, ...)
	} else {
			result = NULL
	}
result
}