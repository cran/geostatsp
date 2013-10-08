

modelRandomFields = function(param){
	

param = fillParam(param)

param["scaleRandomFields"] = param["range"]/2 

if(abs(param["aniso.ratio"]) >  10^(-4)){
	# geometric anisotropy
	scaleTwoDirections = c(1/param["scaleRandomFields"],
			1/(param["aniso.ratio"]*param["scaleRandomFields"]))
	angle = param["aniso.angle.radians"]
	anisoMat = diag(scaleTwoDirections) %*% 
			matrix(c(cos(angle), sin(angle), 
							-sin(angle), cos(angle)),2)
	model=list("$", var=param["variance"],   
			A=anisoMat,
			list("matern", nu=param["rough"]))	
} else {
	model=list("$", var=param["variance"],   
			s=param["scaleRandomFields"],
			list("matern", nu=param["rough"]))	
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
	if(!any(names(param)=="aniso.ratio"))
		param["aniso.ratio"] = 1
	if(!any(names(param)=="aniso.angle.radians")){
		if(any(names(param)=="aniso.angle.degrees")) {
			param["aniso.angle.radians"] = 
					param["aniso.angle.degrees"]*2*pi/360 				
		} else {
			param["aniso.angle.radians"] = 0
		}
	}
	param
}
