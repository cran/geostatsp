
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

    if(length(grep("SpatRaster", class(covariates)))) {
      covariates = covariates[[inModel]]
    } else {
      covariates = covariates[inModel]
    }

    dataFactors = intersect(Sfactor, names(data))

    notInData = setdiff(names(covariates), names(data))

    
    rmethod = rep("bilinear", length(names(covariates)))
    names(rmethod) = names(covariates)
    rmethod[covFactors] = "near"

    
    notLogOffset = ! names(covariates) %in% names(offsetToLogOrig)

    if(any(notLogOffset)){
      if(length(grep("SpatRaster", class(covariates)))) {

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

      covariatesStack = c(cellsSmall, covariatesStack)

      if(length(covariatesForStackData)) {
        covData = stackRasterList(
          covariatesForStackData, 
          data, method=rmethod)
      } else {
        covData = NULL
      }


    } else { # only log offset
      covariatesStack = cellsSmall
      covData = NULL
    }

  for(D in names(offsetToLogOrig)) {
      # loop through offsets which should
      # be aggregated before taking logs
    offsetToLog = covariates[[D]]

    toCrop = union(
      project(
        ext(covariatesStack),
        crs(covariatesStack), 
        crs(offsetToLog)
        ),
      project(ext(data), crs(data),
        crs(offsetToLog)
        )
      )


    offsetToLogCrop = crop(
      offsetToLog, 
      toCrop
      )

    offsetToLogCrop = project(
      offsetToLogCrop,
      y=crs(covariatesStack),
      method='near')

      # aggregate for covariates
    toAggregate = floor(min(res(covariatesStack)/res(offsetToLogCrop)))

    if(any(toAggregate > 1)){
      offsetToLogAgg = aggregate(offsetToLogCrop, fact=toAggregate, 
        fun=sum, na.rm=TRUE)
    } else {
      toAggregate = 1
      offsetToLogAgg = offsetToLogCrop
    }
    offsetToLogAgg = project(offsetToLogAgg, covariatesStack)

    offsetToLogAgg = classify(
      offsetToLogAgg, 
      t(c(-Inf,0,NA)) 
      )

    offsetToLogLogged = log(offsetToLogAgg) - 
    sum(log(rep_len(toAggregate,2)))
    names(offsetToLogLogged) = paste('log',D,sep='')
    covariatesStack = c(covariatesStack, offsetToLogLogged)
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
          gsub("[+]$", "", paste(".~.", 
            paste("offset(log", D, ")", sep=''),
            offsetNotLogged, 
            sep = '+'))
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
      offsetToLogAgg = project(offsetToLogAgg, covariatesStack)
      offsetToLogAgg = classify(
        offsetToLogAgg, 
        t(c(-Inf,0,NA)) 
        )

      offsetToLogLogged = log(offsetToLogAgg) + 
      sum(log(res(covariatesStack))) -
      sum(log(res(offsetToLogCrop)))
      names(offsetToLogLogged) = paste('log',D,sep='')
    }

    covData = c(
      covData,
      offsetToLogLogged
      )

  } # end D in names(offsetToLogOrig)

  covariatesSP = as.points(covariatesStack)
  covariatesDF = values(covariatesSP)

  data = c(data, covData)			


  } else { # if length(inmodel)
    covariatesDF = data.frame()
  }



if(any(res(data)>1.25*res(cellsSmall)))
  warning("data is coarser than grid")

data = c(data, resample(cellsSmall, data, method='near'))	


dataSP = suppressWarnings(as.points(data))
dataDF = values(dataSP)

# get rid of rows with missing response if lgcp with count response

if(names(dataDF)[1] == 'count')
  dataDF = dataDF[!is.na(dataDF$count), ]
  
  # redo factors
# loop through spatial covariates which are factors
for(D in intersect(Sfactor, names(covariatesDF))) {
  theTable = sort(table(dataDF[[D]]), decreasing=TRUE)
  theLevels = levels(covariates[[D]])[[1]]
  if(identical(theLevels, "")) {
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
  covariates, 
  buffer=0) {

  if(missing(covariates)) covariates = list()




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

  # check for missing CRS
  
  if(!nchar(crs(data)) | !nchar(crs(grid)) ) {
    if(nchar(crs(data))) {
      warning("assigning crs of grid to data")
      crs(grid) = crs(data)
    } else if(nchar(crs(grid))) {
      crs(grid) = crs(data)
      warning("assigning crs of data to grid")
    } else if(length(covariates)) {
        if(nchar(crs(covariates[[1]]))) {
          warning("assigning crs of first covariate to data")
          crs(grid) = crs(data) = crs(covariates[[1]])
        } else {
          warning("no crs supplied")
        }
    } else {
      warning("no crs supplied")
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
      names(data)[unlist(lapply(values(data), is.factor))],
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
      rmethod[covFactors] = "near"
      rmethod[intersect(names(covariates), termsInF)] = "near"

      covariatesStack = stackRasterList(
        covariates, 
        template=cellsSmall, 
        method=rmethod)
      covariatesStack = c(cellsSmall, covariatesStack)
      covariatesSP = suppressWarnings(as.points(covariatesStack))
      covariatesDF = values(covariatesSP)
    } else {
      covariatesDF = data.frame()
    }

      # loop through covariates which aren't in data, extract it from `covariates`
    for(D in setdiff(alltermsPlain, names(data))){
      if(is.null(covariates[[D]]))
        warning("cant find covariate '", D, "' in covariates or data")

      if(any(class(covariates[[D]]) == 'SpatRaster')) {
        extractHere = terra::extract(covariates[[D]], 
                                     project(data, crs(covariates[[D]])), ID=FALSE)
      } else {
        extractHere = terra::extract(covariates[[D]], 
                                     project(data, crs(covariates[[D]])))
        extractHere = extractHere[match(1:length(data), extractHere[,'id.y']), 2]
      }
      
      
      if(is.data.frame(extractHere)) {
        if(nrow(extractHere) != nrow(data)) {warning("mismatch in extracted covariates and data")}
        extractHere = extractHere[,grep("^ID$|^id.y", names(extractHere), invert=TRUE), drop=FALSE]
        data[[D]] = extractHere[,1]
      } else {
        data[[D]] = extractHere
      }
    }

    
    # reproject data to grid
    if(!identical(crs(cellsSmall), crs(data))) {
        data = project(data, crs(cellsSmall))
    }
    data$space = suppressWarnings(terra::extract(cellsSmall, data, ID=FALSE, mat=FALSE, dataframe=FALSE))


  # loop through spatial covariates which are factors
    for(D in intersect(Sfactors, names(covariatesDF))) {
      theLevels = levels(covariates[[D]])[[1]]
      idCol = grep("^id$", names(theLevels), ignore.case=TRUE, value=TRUE)[1]
      if(is.na(idCol)) idCol = 1
      if(D %in% names(theLevels)) {
          labelCol = D
      } else {
          labelCol = grep("^category$|^label$", 
            names(theLevels), ignore.case=TRUE, value=TRUE)[1]
      }
      if(is.na(labelCol)) labelCol = 2

      dataD = unlist(data[[D]])
      if(is.factor(dataD)) {

        # give covariatesDF factor levels from data
        if(all(levels(dataD) %in% theLevels[,labelCol])) {
          # match factor levels in data to 
          # factor levels in raster
        theTable = table(dataD)
        # make baseline category the most populous category
        # unless the variable was supplied as a factor in the 'data' argument
        if(D %in% covFactors) {
          theTable = sort(theTable, decreasing=TRUE)
        }
        if(theTable[1]==0) warning("no data in baseline level ", D)
        levelsHave = names(theTable)[theTable > 0]
          covariatesDF[[D]] = factor(
            as.character(covariatesDF[[D]]),
            levels = levelsHave
            )
          data[[D]] = factor(as.character(dataD), levels = levelsHave)
        } else { 
          # levels in data can't be found in raster levels
          # ignore raster levels
          covariatesDF[[D]] = factor(
            covariatesDF[[D]],
            levels = 1:nlevels(data[[D]]),
            labels = levels(data[[D]])
            )

        }

      } else { # data[[D]] isn't a factor

        # choose baseline category 


      theTable = sort(table(data[[D]]), decreasing=TRUE)
      theTable = theTable[theTable > 0]
      if(is.null(theLevels)) {
        theLabels = paste("l", names(theTable),sep="")
      } else {

        if(all(names(theTable) %in% theLevels[,labelCol])) {
      # convert table names to numeric
      # code must work for data where data[[D]] is numeric
          names(theTable) = theLevels[
          match(names(theTable), theLevels[,labelCol])
          , idCol]
        }


        theLabels = as.character(theLevels[
          match(names(theTable), as.character(theLevels[,idCol])), 
          labelCol
          ])
        if(is.numeric(data[[D]])) {
          levelsD = as.integer(names(theTable))
        } else {
          levelsD = theLabels
        }
        if(any(is.na(theLabels))) {
          warning(
            'missing labels in covariate raster ', 
            D, ' level ',
            names(theTable)[is.na(theLabels)][1])
          theLabels[is.na(theLabels)] = 
          names(theTable)[is.na(theLabels)]
        }
      } # end else (not is null thelevels)

      # re-factor data with new baseline category
        data[[D]] = factor(
          data[[D]], 
          levels=levelsD,
          labels=theLabels)     
       covariatesDF[[D]] = factor(
          as.integer(covariatesDF[[D]]), 
          levels=as.integer(names(theTable)),
          labels=theLabels)     
   } # end refactor
   } # end loop D trhoguh factors


   list(
    data=data,
    grid = cellsSmall,
    covariates=covariatesDF
    )
 }
