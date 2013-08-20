
trend.spatial = function(trend, geodata) {
	if(class(geodata) == "SpatialPointsDataFrame")
		geodata = geodata@data
	
	result = model.matrix(trend, data=geodata)
	
	class(result) <- "trend.spatial"
	result
}

