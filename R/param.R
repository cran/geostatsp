

modelRandomFields = function(param){
	

param = fillParam(param)

param["scaleRandomFields"] = param["range"]/2 

if(abs(param["anisoRatio"]) >  10^(-4)){
	# geometric anisotropy
	scaleTwoDirections = c(1/param["scaleRandomFields"],
			1/(param["anisoRatio"]*param["scaleRandomFields"]))
	angle = param["anisoAngleRadians"]
	anisoMat = diag(scaleTwoDirections) %*% 
			matrix(c(cos(angle), sin(angle), 
							-sin(angle), cos(angle)),2)
	model=list("$", var=param["variance"],   
			A=anisoMat,
			list("matern", nu=param["shape"]))	
} else {
	model=list("$", var=param["variance"],   
			s=param["scaleRandomFields"],
			list("matern", nu=param["shape"]))	
}	

model 
}

fillParam = function(param) {
	names(param) = gsub("^var$", "variance", names(param))
	
	if(!any(names(param)=="variance") & any(names(param)=="sdSpatial"))
		param["variance"]= param["sdSpatial"]^2
	
# if still no variance set it to 1.
	if(!any(names(param)=="variance")) {
		param["variance"] = 1
	}
	
# fill in anisotropy parameters	
	if(!any(names(param)=="anisoRatio"))
		param["anisoRatio"] = 1
	if(!any(names(param)=="anisoAngleRadians")){
		if(any(names(param)=="anisoAngleDegrees")) {
			param["anisoAngleRadians"] = 
					param["anisoAngleDegrees"]*2*pi/360 				
		} else {
			param["anisoAngleRadians"] = 0
		}
	}
	param
}
