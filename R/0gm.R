
allVarsP = function(formula) {
  # return vector of variable names in the formula
  allterms = rownames(attributes(terms(formula))$factors)
  if(!length(allterms)) {
  	# might be intercepty only but there may be an offset
  	# add a junk variable so that the terms function produces a factors
  	allterms =  setdiff(
  	rownames(attributes(terms(
  		update.formula(formula, .~ .+ junk)
  		))$factors),
  		'junk')
  }
  
  firstTerm = as.character(formula)
  firstTerm = trimws(firstTerm[-c(1, length(firstTerm))])
  allterms = setdiff(allterms, firstTerm)
  
  # look for values=1:stuff in inla formula and replace with seq
  
  if(length(allterms)) {
    inlaValuesPattern = 'values[[:space:]]+?=[[:space:]]+?'
    haveInlaValues = grep(inlaValuesPattern, allterms, value=TRUE)
    notInlaValues = grep(inlaValuesPattern, allterms, 
      value=TRUE, invert=TRUE)
    allterms = c(
      unique(unlist(strsplit(notInlaValues, ":"))),
      haveInlaValues
    )
  }
  allterms = gsub("[[:space:]]", "", allterms)
# remove offset( or factor(
  alltermsPlain = gsub("^[[:alpha:]]+\\(|\\)$|[,].*", "", allterms)
  attributes(alltermsPlain)$orig = allterms
  alltermsPlain
}

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
  
  allterms = allVarsP(formula)
  
  alltermsFull = attributes(allterms)$orig

  # find factors
  
  theFactors = grep("^factor", alltermsFull, value=T)
  theFactors = gsub("^factor\\(|\\)$", "", theFactors)
  
  termsInF = grep("^f[(]", alltermsFull, value=TRUE)
  termsInF = gsub("^f[(]|[,].*|[[:space:]]", "", termsInF)
  
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
  
  inModel = gsub("^[[:alnum:]]+[(]|[,].*|[[:space:]]",
    "",  alltermsFull)
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
      '', offsetToLogOrig
    )
  }
  
  Sfactor = c(
    dataFactors,
    covFactors,
    theFactors
  )
  covFactors = intersect(Sfactor,names(covariates))
  dataFactors = intersect(Sfactor,names(data))
  
  inModel = intersect(inModel, names(covariates))
  if(length(inModel)) {	
    if(length(grep("^Raster", class(covariates)))) {
      covariates = covariates[[inModel]]
    } else {
      covariates = covariates[inModel]
    }
    dataFactors = intersect(Sfactor, names(data))
    notInData = setdiff(names(covariates), names(data))
    
    rmethod = rep("bilinear", length(names(covariates)))
    names(rmethod) = names(covariates)
    rmethod[covFactors] = "ngb"
    
    notLogOffset = ! names(covariates) %in% names(offsetToLogOrig)
    if(any(notLogOffset)){
      if(length(grep("^Raster", class(covariates)))) {
        covariatesForStack = covariates[[which(notLogOffset)]]
        covariatesForStackData = 
          covariates[[notInData]]
      } else {
        covariatesForStack = covariates[notLogOffset]
        covariatesForStackData = 
          covariates[notInData]
      }
      covariatesStack = stackRasterList(
        covariatesForStack,
        cellsSmall, method=rmethod)
      
      covariatesStack = stack(cellsSmall, covariatesStack)
      covData = stackRasterList(
        covariatesForStackData, 
        data, method=rmethod)
      
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
        crs=crs(covariatesStack),
        method='ngb')
      
      # aggregate for covariates
      toAggregate = floor(min(res(covariatesStack)/res(offsetToLogCrop)))
      
      if(any(toAggregate > 1)){
        offsetToLogAgg = aggregate(offsetToLogCrop, fact=toAggregate, 
          fun=sum, na.rm=TRUE)
      } else {
        toAggregate = 1
        offsetToLogAgg = offsetToLogCrop
      }
      offsetToLogAgg = projectRaster(offsetToLogAgg, covariatesStack)
      
      offsetToLogAgg = reclassify(
        offsetToLogAgg, 
        t(c(-Inf,0,NA)) 
      )
      
      offsetToLogLogged = log(offsetToLogAgg) - 
        sum(log(rep_len(toAggregate,2)))
      names(offsetToLogLogged) = paste('log',D,sep='')
      covariatesStack = stack(covariatesStack, offsetToLogLogged)
      toDrop = which(alltermsFull==offsetToLogOrig[D])
      
      # the offsets
      allOffsets = grep(
        "^offset\\([[:print:]]+\\)$", 
        gsub("[[:space:]]+", "", alltermsFull), value=TRUE)
      offsetNotLogged = grep(
        "^offset\\([[:print:]]+,log=TRUE\\)$", 
        gsub("[[:space:]]+", "", allOffsets), 
        invert=TRUE, value=TRUE)
      
      
      # any other offsets would also have been removed
      # by drop.terms
      
      
      if(length(allOffsets)< length(alltermsFull)) {
        formula = update.formula(
          drop.terms(terms(formula), dropx=toDrop, keep.response=TRUE),
          as.formula(
            paste(".~.", 
              paste("offset(log", D, ")", sep=''),
              offsetNotLogged, 
              sep = '+')
          ) 	
        )
      } else {
        # only offsets in the model
        # drop.terms won't work
        formula = update.formula(
          formula,
          as.formula(
            paste(".~", 
              paste("offset(log", D, ")", sep=''),
              offsetNotLogged, 
              sep = '+')
          )
        
        )	
      }
      
      rmethod[paste('log',D,sep='')] = 'bilinear'
      
      # aggregate for data
      toAggregateData = floor(min(res(data)/res(offsetToLogCrop)))
      if(toAggregateData != toAggregate & toAggregateData > 1 ){
        offsetToLogAgg = aggregate(offsetToLogCrop, fact=toAggregateData, fun=sum)
        offsetToLogAgg = projectRaster(offsetToLogAgg, covariatesStack)
        offsetToLogAgg = reclassify(
          offsetToLogAgg, 
          t(c(-Inf,0,NA)) 
        )
        
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


#############
# data is a SpatialPointsDataFrame
#############

gm.dataSpatial = function(
  formula, data,  grid, 
  covariates=NULL, 
  buffer=0) {
  
  # check response variable is in data
  if(!all.vars(formula)[1] %in% names(data)){
    warning(paste(
        'response variable',
        all.vars(formula)[1],
        'not found in data'
      ))
  }
  
  alltermsPlain = allVarsP(formula)
  
# find factors
  allterms = attributes(alltermsPlain)$orig
  
  # remove covariates not in the model
  keepCovariates = intersect(alltermsPlain, names(covariates))
  if(is.list(covariates)) {
    covariates = covariates[keepCovariates]
  } else {
    if(length(keepCovariates)) {
      covariates = covariates[[keepCovariates]]	
    } else {
      covariates = list()
    }
  }
  
  theFactors = grep("^factor", allterms, value=T)
  theFactors = gsub("^factor\\(|\\)$", "", theFactors)
  
  
  termsInF = grep("^f[(]", allterms, value=TRUE)
  termsInF = gsub("^f[(]|[,].*|[[:space:]]", "", termsInF)
  
  
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
    rmethod[intersect(names(covariates), termsInF)] = "ngb"
    
    
    covariatesStack = stackRasterList(
      covariates, 
      template=cellsSmall, 
      method=rmethod)
    covariatesStack = stack(cellsSmall, covariatesStack)
    
    
    covariatesSP = as(covariatesStack, "SpatialPointsDataFrame")
    covariatesDF = covariatesSP@data
  } else {
    covariatesDF = data.frame()
  }
  
  # loop through covariates which aren't in data, extract it from `covariates`
  for(D in setdiff(alltermsPlain, names(data))){
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
  data$space = suppressWarnings(extract(cellsSmall, data))
  
  # loop through spatial covariates which are factors
  for(D in intersect(Sfactors, names(covariatesDF))) {
    theTable = sort(table(data[[D]]), decreasing=TRUE)
    theLevels = levels(covariates[[D]])[[1]]
    if(is.null(theLevels)) {
      theLabels = paste("l", names(theTable),sep="")
    } else {
      idCol = grep("^id$", names(theLevels), ignore.case=TRUE)[1]
      if(is.na(idCol)) idCol = 1
      labelCol = grep("^category$|^label$", names(theLevels), ignore.case=TRUE)[1]
      if(is.na(labelCol)) labelCol = 2
 
    if(all(names(theTable) %in% theLevels[,labelCol])) {
      # convert table names to numeric
      # code must work for data where data[[D]] is numeric
      names(theTable) = theLevels[
        match(names(theTable), theLevels[,labelCol])
        , idCol]
    }

      theLabels = as.character(theLevels[match(names(theTable), theLevels[,idCol]), labelCol])
      if(any(is.na(theLabels))) {
        warning(
          'missing labels in covariate raster ', 
          D, ' level ',
          names(theTable)[is.na(theLabels)][1])
        theLabels[is.na(theLabels)] = 
          names(theTable)[is.na(theLabels)]
      }
    }
    data[[D]] = factor(
      as.integer(data[[D]]), 
      levels=as.integer(names(theTable)),
      labels=theLabels)			
    covariatesDF[[D]] = factor(
      as.integer(covariatesDF[[D]]), 
      levels=as.integer(names(theTable)),
      labels=theLabels)			
  }
  
  
  list(
    data=data,
    grid = cellsSmall,
    covariates=covariatesDF
  )
}
