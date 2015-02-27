
#+ brickFunction
spdfToBrick = function(x, 
    template,
    logSumExpected=FALSE,
    pattern = '^logExpected_surfaceArea_[[:digit:]]+$'
) {
  
  if(logSumExpected){
    pattern = gsub("\\^logE", "e", pattern)
  }
  
  
  
  if(class(x)=='SpatialPolygonsDataFrame'){
    x = list(x) 
  }
  if(is.null(names(x)))
    names(x) = as.character(1:length(x))
  
  if(is.numeric(template)) {
    template = squareRaster(
        x[[1]], template[1])
  }
  
  forRaster = NULL
  
  haveRgdal = requireNamespace('rgdal', quietly=TRUE)
  
  for(Dcensus in rev(names(x))){
    
    if(haveRgdal &
        !identical(
            projection(x[[Dcensus]]),
            projection(template[[Dcensus]])
            )){
      x[[Dcensus]] = spTransform(
          x[[Dcensus]],
          CRS(projection(template))
      )
    }
    
    Sid = rasterize(
        x[[Dcensus]],
        template,
        field=1:length(x[[Dcensus]])
    )
    
    forRaster = cbind(
        as.matrix(x[[Dcensus]]@data[,
                grep(pattern, names(x[[Dcensus]]))
            ])[values(Sid),
        ],
        forRaster)
  }
  
  if(logSumExpected){
    forRaster = apply(forRaster, 1, sum, na.rm=TRUE)
    forRaster[forRaster<=0] = NA
    forRaster = matrix(log(forRaster))
  }
  
  if(is.matrix(forRaster)){
    result = brick(template, nl=ncol(forRaster))
    values(result) = forRaster
  } else if(is.vector(forRaster)){
    result = template
    values(result) = forRaster
  } else {
    warning("no data extracted")
    result = template
  }  
  
  result
}

#'

