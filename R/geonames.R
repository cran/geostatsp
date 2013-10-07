GNcities = function(north, east, south, west, lang = "en", maxRows = 10) {
	
	
	if(is.vector(north) ) {
		theproj = NULL
		
	} else {
		
		theproj = try(proj4string(north),silent=TRUE)
		if(class(theproj)=="try-error")
			theproj = NULL

		# do this because bbox(mybbox) != mybbox
		# but bbox(extent(mybbox) = mybbox
		north = bbox(extent(north))

		
		if(!is.null(theproj)) {
			# transform to long-lat
			north = SpatialPoints(t(north),
					proj4string=CRS(theproj))
			
			north = bbox(spTransform(
							north, CRS("+proj=longlat")
					))			
		}
		
		east = north[1,2]
		west = north[1,1]
		south = north[2,1]
		north=north[2,2]
		
	}

	result = geonames::GNcities(north,east,south,west,lang,maxRows)
	result = SpatialPointsDataFrame(result[,c("lng","lat")], data=result, 
			proj4string=CRS("+proj=longlat"))

	if(!is.null(theproj))
		result = spTransform(result, CRS(theproj))
	
	result
}