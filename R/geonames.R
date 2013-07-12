GNcities = function(north, ...) {
	UseMethod("GNcities")
	
}

GNcities.matrix = function(north, ...){
	if(!all(dim(north)==2)){
		warning("wrong dimensions, is north a bounding box?")
	}
	result = 	geonames::GNcities(north[2,2], north[1,2], north[2,1], north[1,1], ...)
	SpatialPointsDataFrame(result[,c("lng","lat")], data=result, 
			proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
}

GNcities.Extent = function(north, ...){
	result=geonames::GNcities(
			north@ymax, north@xmax, 
			north@ymin, north@xmin, 
			...)
	SpatialPointsDataFrame(result[,c("lng","lat")], data=result, 
			proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
	
}

GNcities.default = function(north, ...){
	result=geonames::GNcities(north, ...)			
	SpatialPointsDataFrame(result[,c("lng","lat")], data=result, 
			proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
	
}