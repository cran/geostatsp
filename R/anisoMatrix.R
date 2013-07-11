anisoMatrix = function(angle, range=NULL, scale=1/range) {
	
	if(length(scale)!= 2)
		warning("scale should be length 2")
	
	result = diag(scale[1:2]) %*% 
			matrix(c(cos(angle), sin(angle), 
							-sin(angle), cos(angle)),2)
	
	result
	
}