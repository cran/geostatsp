openmap = function(upperLeft, ...) {
	UseMethod("openmap")
	
}

openmap.matrix = function(upperLeft, ...){
	if(!all(dim(upperLeft)==2)){
		warning("wrong dimensions, is xlim a bounding box?")
	}
	OpenStreetMap::openproj(
	OpenStreetMap::openmap(
			upperLeft=c(upperLeft[2,2],upperLeft[1,1]), 
			lowerRight=c(upperLeft[2,1],upperLeft[1,2]),
		 ...)
	)
	 
}

openmap.Extent = function(upperLeft, ...){
	OpenStreetMap::openproj(
			OpenStreetMap::openmap(
			c(upperLeft@ymax, upperLeft@xmin), 
			c(upperLeft@ymin, upperLeft@xmax), 
			...)
)
}

openmap.default = function(upperLeft, ...){
	OpenStreetMap::openproj(
			OpenStreetMap::openmap(upperLeft, ...)			
	)
}
 