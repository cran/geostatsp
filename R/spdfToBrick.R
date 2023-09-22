
#+ brickFunction
spdfToBrick = function(x, 
    template,
    logSumExpected=FALSE,
    pattern = '^expected_[[:digit:]]+$'
) {
  
  
  
  if(any(class(x)=='SpatVector')){
		if(ncol(x)==1) {
			pattern=names(x)
		}
    x = list(x) 
  }
  if(is.null(names(x)))
    names(x) = as.character(1:length(x))
  
  if(is.numeric(template)) {
    template = squareRaster(
        x[[1]], template[1])
  }
  
  forRaster = NULL
  

  for(Dcensus in rev(names(x))){
    
    if(!identical(
            crs(x[[Dcensus]]),
            crs(template)
            )){
      x[[Dcensus]] = project(
          x[[Dcensus]],
          crs(template)
      )
    }
    
		Sx = 1:length(x[[Dcensus]])
		
    Sid = rasterize(
        x[[Dcensus]],
        template,
        field=Sx
    )
		
		# number of raster cells assigned to each polygon
		nCellPerPoly = rep(1, length(Sx))
		nCellTable = table(values(Sid))
		nCellPerPoly[as.numeric(names(nCellTable))] = nCellTable
		
		dataHere = as.matrix(
			values(
				x[[Dcensus]]
			)[,
        grep(pattern, names(x[[Dcensus]])), drop=FALSE
    ]) / nCellPerPoly

		dataHere[is.na(dataHere)] = 0

		# assign expected counts to raster cells
		forRasterHere = dataHere[values(Sid), ,drop=FALSE]
		colnames(forRasterHere) = colnames(dataHere)
		forRasterHere[is.na(forRasterHere)] = 0
		
		# polygons not assigned to cells
		notInRaster = which(! Sx %in% values(Sid))
		if(length(notInRaster)){
			polyCentres = centroids(x[[Dcensus]][notInRaster,])
			polyCell = cellFromXY(template, crds(polyCentres))
			
			dataNotInRaster = aggregate(
					dataHere[notInRaster,,drop=FALSE], 
					list(cell=polyCell), 
					FUN=sum, na.rm=TRUE)
			forRasterHere[dataNotInRaster$cell, ] = 
					forRasterHere[dataNotInRaster$cell, ,drop=FALSE] + 
					as.matrix(dataNotInRaster[,-1,drop=FALSE])
		}
		
    forRaster = cbind(
				forRasterHere,
        forRaster)
  }
  

	
  if(logSumExpected){
    forRaster = apply(forRaster, 1, sum, na.rm=TRUE)
    forRaster[forRaster<=0] = NA	
		# divide by cell size to get intensity
    forRaster = matrix(log(forRaster) - sum(log(res(template))))
		colnames(forRaster) = 'logExpected'
  } else {
		# divide by cell size to get intensity
    forRaster = forRaster / prod(res(template))
	}
  
  if(is.matrix(forRaster)){
    result = rast(template, nlyrs=ncol(forRaster))
    terra::values(result) = forRaster
    names(result) = colnames(forRaster)
  } else if(is.vector(forRaster)){
    result = template
    terra::values(result) = forRaster
    names(result) = 'expected'
  } else {
    warning("no data extracted")
    result = template
  }  
  
  result
}

#'

