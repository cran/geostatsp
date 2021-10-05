lgcp = function(formula=NULL, data,  grid, covariates=NULL, 
	border,
	...) {
	
	if(!missing(border)){
		if(!.compareCRS(data, border))
			border = spTransform(border, data@proj4string)
	}
	
	if(is.numeric(grid)) {
		if(!missing(border)){
			cells = squareRaster(border,grid)
		} else {
			cells = squareRaster(data,grid)
		}
	} else {
		cells = squareRaster(grid)
	}
	
# create data
	
	if(!missing(border)) {
		inBorder = over(
			data, 
			as(border, 'SpatialPolygons')
			)
		data = data[!is.na(inBorder),]
	}
	
	counts = rasterize(
		SpatialPoints(data), 
		cells, fun="count")
	names(counts) = "count"
	counts[is.na(counts)] = 0

	if(!missing(border)) {
		counts = mask(counts, border)
	}
	
	
# the formula	
	if(is.null(formula)) {
		formula = as.formula(
			paste(c("count ~ 1", names(covariates)), collapse="+")
			)
	}

	formula	= update.formula(formula,
		.~.+offset(logCellSize) 
		)
	formula = update.formula(formula, count ~ .)

	alltermsFull = rownames(attributes(terms(formula))$factors)[-1]
	offsetToLogOrig = grep(
		"^offset\\([[:print:]]+,log=TRUE\\)$", 
		gsub("[[:space:]]+", "", alltermsFull))
	offsetToLogOrig = alltermsFull[offsetToLogOrig]

	offsetToMask = gsub(
		"^[[:space:]]?offset\\(|,[[:space:]]?log[[:space:]]?=[[:space:]]?TRUE[[:space:]]?\\)[[:space:]]?$",
		'', offsetToLogOrig
		)[1]

	if(!missing(border)) {
		# set values of the offset to zero outside the border
#	instead of masking as was done formerly with	counts = raster::mask(counts, border)
		
		if(length(offsetToLogOrig)) {
			if(offsetToMask %in% names(covariates)) {
				if(!.compareCRS(covariates[[offsetToMask]], border)) {
					borderM = spTransform(border, 
						covariates[[offsetToMask]]@crs)
				} else {
					borderM = border
				}
				covariates[[offsetToMask]] = raster::mask(
					covariates[[offsetToMask]], borderM
					)
			}
		}
	}
	
	
	# cell size offset

	if(length(grep("^Raster", class(covariates)))) {
		# add a raster layer for log cell size
		logCellSize = raster(covariates)
		values(logCellSize) = sum(log(res(cells)))
		names(logCellSize) = 'logCellSize'
		covariates = addLayer(logCellSize, covariates)
	} else {
		# create a raster and put it in the covariate list
		logCellSize = cells
		names(logCellSize) = "logCellSize"
		values(logCellSize) =  sum(log(res(cells)) )
		if(length(covariates)){
			covariates = c(covariates, logCellSize=logCellSize)
		} else {
			covariates = logCellSize
		}
	}
	
	dots = list(...)
	if(!any(names(dots)=='family')) {
		dots$family='poisson'
	}
	dots$formula = formula
	dots$data = counts
	dots$grid = cells
	dots$covariates=covariates
	

	result = do.call(glgm, dots)

	result

}


