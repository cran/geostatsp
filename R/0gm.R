gm.dataRaster = function(
    formula,
    data, grid=data,
    covariates=NULL,
    buffer=0){

  if(abs(diff(res(grid)))>0.000001 )
    warning("data is not on a square grid")
  
  cellsBoth = cellsBuffer(grid, buffer)			
  cellsSmall = cellsBoth$small
	
	# find factors
	
	alltermsFull = rownames(attributes(terms(formula))$factors)[-1]
  if(!length(alltermsFull)){
    # maybe there's offsets
    # this thing works if that's the case
    alltermsFull = as.character(attributes(terms(formula))$variables)[-(1:2)]
  }
	allterms = grep(":", alltermsFull, invert=TRUE, value=TRUE)
	allterms = gsub("[[:space:]]", "", allterms)
	
	theFactors = grep("^factor", allterms, value=T)
	theFactors = gsub("^factor\\(|\\)$", "", theFactors)
	
	
	covFactors = NULL
	for(D in names(covariates)) {
		if(is.factor(covariates[[D]]))
			covFactors = c(D, covFactors)
	}
	dataFactors = NULL
	for(D in names(data)) {
		if(is.factor(data[[D]]))
			dataFactors = c(D, dataFactors)
	}
	
  inModel = gsub("^[[:alnum:]]+\\(",
      "",
      alltermsFull)
  inModel = gsub(
      "(,([[:alnum:]]|=|[[:space:]])+)?\\)$",#+?[[:space:]]?\\)[[:space:]]?$",
      "", inModel)
  
  offsetToLogOrig = grep(
      "^offset\\([[:print:]]+,log=TRUE\\)$", 
      gsub("[[:space:]]+", "", alltermsFull))
  offsetToLogOrig = alltermsFull[offsetToLogOrig]
  if(length(offsetToLogOrig)) {
  names(offsetToLogOrig) = gsub(
      "^[[:space:]]?offset\\(|,[[:space:]]?log[[:space:]]?=[[:space:]]?TRUE[[:space:]]?\\)[[:space:]]?$",
      '', offsetToLogOrig)
  }
  
	Sfactor = c(
		dataFactors,
		covFactors,
		theFactors
	)
	covFactors = intersect(Sfactor,names(covariates))
	dataFactors = intersect(Sfactor,names(data))
	
  inModel = intersect(inModel, names(covariates))
  covariates = covariates[inModel]
	if(length(inModel)) {		
    dataFactors = intersect(Sfactor, names(data))
    notInData = setdiff(names(covariates), names(data))
    
		rmethod = rep("bilinear", length(names(covariates)))
		names(rmethod) = names(covariates)
		rmethod[covFactors] = "ngb"
		
    notLogOffset = ! names(covariates) %in% names(offsetToLogOrig)
    if(any(notLogOffset)){
		covariatesStack = stackRasterList(
        covariates[notLogOffset],
        cellsSmall, method=rmethod)

    covariatesStack = stack(cellsSmall, covariatesStack)
    covData = stackRasterList(covariates[notInData], data, method=rmethod)
    
  } else {
    covariatesStack = cellsSmall
    covData = NULL
  }
  
    
    for(D in names(offsetToLogOrig)) {
      # loop through offsets which should
      # be aggregated before taking logs
      offsetToLog = covariates[[D]]
      
      toCrop = merge(
          projectExtent(covariatesStack, 
            crs(offsetToLog)
    ),
          projectExtent(data, 
              crs(offsetToLog)
          )
      )
      
      
      offsetToLogCrop = raster::crop(
          offsetToLog, 
          toCrop
      )
      
      offsetToLogCrop = projectRaster(
          offsetToLogCrop,
          crs=crs(covariatesStack))

      # aggregate for covariates
      toAggregate = floor(min(res(covariatesStack)/res(offsetToLogCrop)))
      if(toAggregate > 1){
        offsetToLogAgg = aggregate(offsetToLogCrop, fact=toAggregate, fun=sum)
      }
      offsetToLogAgg = projectRaster(offsetToLogAgg, covariatesStack)
      offsetToLogLogged = log(offsetToLogAgg) + 
          sum(log(res(covariatesStack))) -
          sum(log(res(offsetToLogCrop)))
      names(offsetToLogLogged) = paste('log',D,sep='')
      covariatesStack = stack(covariatesStack, offsetToLogLogged)
      toDrop = which(alltermsFull==offsetToLogOrig[D])
#      drop.terms(terms(formula), toDrop, keep.response=TRUE)
    # doesn't work, will drop other offsets

    formula = as.formula(gsub(offsetToLogOrig[D],
        paste("offset(log",D,")",sep=''),
        format(formula),
        fixed=TRUE))
    
    rmethod[paste('log',D,sep='')] = 'bilinear'

    # aggregate for data
toAggregateData = floor(min(res(data)/res(offsetToLogCrop)))
if(toAggregateData != toAggregate & toAggregateData > 1 ){
  offsetToLogAgg = aggregate(offsetToLogCrop, fact=toAggregateData, fun=sum)
  offsetToLogAgg = projectRaster(offsetToLogAgg, covariatesStack)
  offsetToLogLogged = log(offsetToLogAgg) + 
      sum(log(res(covariatesStack))) -
      sum(log(res(offsetToLogCrop)))
  names(offsetToLogLogged) = paste('log',D,sep='')
}


    covData = stack(
        covData,
        offsetToLogLogged
        )
    
    } # end in names(offsetToLogOrig)
    covariatesSP = as(covariatesStack, "SpatialPointsDataFrame")
    covariatesDF = covariatesSP@data
    
    data = stack(data, covData)			
    
    
	} else {
		covariatesDF = data.frame()
	}
	
	if(any(res(data)>1.25*res(cellsSmall)))
		warning("data is coarser than grid")
	
	data = stack(data, resample(cellsSmall, data, method='ngb'))	
	dataSP = as(data, "SpatialPointsDataFrame")
	dataDF =dataSP@data

	# redo factors
# loop through spatial covariates which are factors
for(D in intersect(Sfactor, names(covariatesDF))) {
	theTable = sort(table(dataDF[[D]]), decreasing=TRUE)
	theLevels = levels(covariates[[D]])[[1]]
	if(is.null(theLevels)) {
		theLabels = paste("l", names(theTable),sep="")
	} else {
		theLabels = theLevels[
				match(as.integer(names(theTable)), theLevels$ID)
				,"Category"]
	}
	dataDF[[D]] = factor(dataDF[[D]], levels=as.integer(names(theTable)),
			labels=theLabels)			
	covariatesDF[[D]] = factor(covariatesDF[[D]], levels=as.integer(names(theTable)),
			labels=theLabels)			
	
}

  list(
    data=dataDF,
    grid=cellsSmall,
    covariates=covariatesDF,
    formula = formula
    )
}

gm.dataSpatial = function(
    formula, data,  grid, 
		covariates=NULL, 
		buffer=0) {

# find factors
	allterms = colnames(attributes(terms(formula))$factors)
	allterms = grep(":", allterms, invert=TRUE, value=TRUE)
	allterms = gsub("[[:space:]]", "", allterms)


	
	theFactors = grep("^factor", allterms, value=T)
	theFactors = gsub("^factor\\(|\\)$", "", theFactors)

	covFactors = NULL
	for(D in names(covariates)) {
		if(is.factor(covariates[[D]]))
			covFactors = c(D, covFactors)
	}

	
	Sfactors = c(
			names(data)[unlist(lapply(data@data, is.factor))],
			covFactors,
			theFactors
	)
	Sfactors = unique(Sfactors)
	covFactors = intersect(Sfactors,names(covariates))
	
	cantFind = setdiff(Sfactors, c(names(data), names(covariates)))
	if(length(cantFind))
		warning("can't find variables", cantFind)
	
	# the grid
	cellsBoth = cellsBuffer(grid, buffer)			
	cellsSmall = cellsBoth$small
  
	# 
	if(length(names(covariates))) {
 
		dataFactors = intersect(Sfactors, names(data))
		
		rmethod = rep("bilinear", length(names(covariates)))
		names(rmethod) = names(covariates)
		rmethod[covFactors] = "ngb"
		
		
		covariatesStack = stackRasterList(covariates, template=cellsSmall, method=rmethod)
		covariatesStack = stack(cellsSmall, covariatesStack)
 
		
		covariatesSP = as(covariatesStack, "SpatialPointsDataFrame")
		covariatesDF = covariatesSP@data
	} else {
		covariatesDF = data.frame()
	}
	
	# loop through factors which aren't in data, extract it from covariates
	for(D in setdiff(all.vars(formula), names(data))){
		if(is.null(covariates[[D]]))
			warning("cant find covariate '", D, "' in covariates or data")
		if(!.compareCRS(covariates[[D]], data, unknown=TRUE) ) {
			if(requireNamespace('rgdal', quietly=TRUE ) ) { 
				data[[D]] = raster::extract(covariates[[D]], 
						spTransform(data, CRSobj=CRS(projection(covariates[[D]]))))
			} else warning("need rgdal if covariates and data are different projections")
		} else {
			data[[D]] = raster::extract(covariates[[D]], 
					data) 
		}
	}
	data$space = extract(cellsSmall, data) 

	# loop through spatial covariates which are factors
	for(D in intersect(Sfactors, names(covariatesDF))) {
		theTable = sort(table(data[[D]]), decreasing=TRUE)
		theLevels = levels(covariates[[D]])[[1]]
		if(is.null(theLevels)) {
			theLabels = paste("l", names(theTable),sep="")
		} else {
			idCol = grep("^id$", names(theLevels), ignore.case=TRUE)[1]
			if(!length(idCol)) idCol = 1
			labelCol = grep("^category$|^label$", names(theLevels), ignore.case=TRUE)[1]
			if(!length(labelCol)) labelCol = 2
			
			theLabels = theLevels[
				match(as.integer(names(theTable)), theLevels[,idCol])
				,labelCol]
		}
		data[[D]] = factor(data[[D]], levels=as.integer(names(theTable)),
				labels=theLabels)			
		covariatesDF[[D]] = factor(covariatesDF[[D]], levels=as.integer(names(theTable)),
				labels=theLabels)			
		
	}


list(
		data=data,
		grid = cellsSmall,
		covariates=covariatesDF
	)
}
