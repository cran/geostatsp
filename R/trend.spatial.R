
trend.spatial = function(trend, geodata) {
	if(class(geodata) == "SpatialPointsDataFrame")
		geodata = geodata@data
	
	result = model.matrix(trend, data=geodata)
	
	class(result) <- "trend.spatial"
	result
}


"trend.spatialOld" <-
		function (trend, geodata, add.to.trend) 
{  
#	if(!missing(geodata)){
	#	if(any(class(geodata) %in% ls(pattern=glob2rx("Spatial*DataFrame"), pos="package:sp")))
	#		geodata <- geodata@data
#		attach(geodata, pos=2)
#		if(!is.null(geodata$covariate)){
	#		attach(geodata$covariate, pos=3)
	#		on.exit(detach("geodata$covariate"), add=TRUE)
#		} 	
#		on.exit(detach("geodata"), add=TRUE)
#	}
	if(inherits(trend, "formula")) {
		#    require(methods)
		#    if(exists("trySilent")){
		if(!missing(geodata)){
			trend.mat <- try(
					model.matrix(trend, data=geodata),
					silent=TRUE) 
	} else {
		trend.mat <- try(model.matrix(trend),silent=TRUE)
	}
		#    }
		#    else{
		#      error.now <- options()$show.error.messages
		#      if (is.null(error.now) | error.now) 
		##      options(show.error.messages = FALSE)
		#     trend.mat <- try(model.matrix(trend))
		#   }    
		if (inherits(trend.mat, "try-error")) 
			stop("\ntrend elements not found")
	}
	else {
		if(class(geodata) == "data.frame") {
			Npoints = dim(geodata)[1]
		} else if(class(geodata) == "geodata") {
			Npoints = dim(geodata$coords)[1]
		} else if (is.vector(geodata)) {
			Npoints = geodata[1]
		} else Npoints = NULL
		
		if(mode(trend) == "numeric")
			trend.mat <- unclass(trend)
		else if (trend == "zeros"){
			if(missing(geodata))
				stop("argument geodata must be provided with trend=\"cte\"")
			trend.mat <- as.matrix(rep(0, Npoints))
		}
		else if (trend == "cte"){
			if(missing(geodata))
				stop("argument geodata must be provided with trend=\"cte\"")
			trend.mat <- as.matrix(rep(1, Npoints))
		}
		else if (trend == "1st"){
			if(missing(geodata))
				stop("argument geodata must be provided with trend=\"1st\"")
			trend.mat <- cbind(1, geodata$coords)
		}
		else if (trend == "2nd"){ 
			if(missing(geodata))
				stop("argument geodata must be provided with trend=\"2nd\"")
			trend.mat <- cbind(1, geodata$coords, geodata$coords[,1]^2,
					geodata$coords[,2]^2,
					geodata$coords[,1] * geodata$coords[,2])
		}
		else stop("external trend must be provided for data locations to be estimated using the arguments trend.d and trend.l. Allowed values are the strings \"cte\", \"1st\", \"2nd\" or  a model formula")
	}

	trend.mat <- as.matrix(trend.mat)
	if(!missing(add.to.trend)){
		if(missing(geodata))
			trend.mat <- cbind(trend.mat, trend.spatial(add.to.trend)[,-1])
		else
			trend.mat <- cbind(trend.mat, trend.spatial(add.to.trend, geodata = geodata)[,-1])
	}
#	dimnames(trend.mat) <- list(NULL, NULL)
	oldClass(trend.mat) <- "trend.spatial"
	return(trend.mat)
}
