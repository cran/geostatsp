
gm.nullFormula = 
		function(formula=NULL, data, grid, 
		covariates=NULL, ...) {
	formula =  1 
	callGeneric(
			formula=formula, data=data,
			grid=grid,
			covariates=covariates, 
			...
	)	
}

gm.numericFormula = 
		function(formula, data, grid, 
				covariates=NULL, ...) {
	
formula = names(data)[formula]
callGeneric(
		formula=formula, data=data,
		grid=grid,
		covariates=covariates, 
		...
)	
}

gm.characterFormula = 
		function(formula, data, grid, 
				covariates=NULL, ...) {


	if(length(names(covariates)))
	names(covariates) = gsub("[[:punct:]]|[[:space:]]","_", names(covariates))
if(length(covariates) & !length(names(covariates))) 
	names(covariates) = paste("c", 1:length(covariates),sep="")			

if(length(formula)==1)
	formula = unique(c(formula, names(covariates)))

formula = paste(formula[1] , "~",
		paste(formula[-1], collapse=" + ")
)
formulaAsFormula = as.formula(formula)

callGeneric(
		formula=formulaAsFormula, data=data,
		grid=grid,
		covariates=covariates, 
		...)
}		


gm.gridNumeric = function(formula, data, grid, covariates=NULL, ...) {
		

	gridRaster = squareRaster(data, grid)
	
	callGeneric(
			formula=formula, data=data,
			grid = gridRaster,
			covariates=covariates, 
			...
	)	
}

gm.dataRaster = function(formula,data,grid, covariates=NULL, buffer=0,...){

	
	cellsBoth = cellsBuffer(grid, buffer)			
	cellsSmall = cellsBoth$small
	
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
	dataFactors = NULL
	for(D in names(data)) {
		if(is.factor(data[[D]]))
			dataFactors = c(D, dataFactors)
	}
	
	Sfactor = c(
		dataFactors,
		covFactors,
		theFactors
	)
	covFactors = intersect(Sfactor,names(covariates))
	dataFactors = intersect(Sfactor,names(data))
	
	if(length(covariates)) {
		dataFactors = intersect(Sfactor, names(data))
		
		rmethod = rep("bilinear", length(names(covariates)))
		names(rmethod) = names(covariates)
		rmethod[covFactors] = "ngb"
		
		covariatesStack = stackRasterList(covariates, cellsSmall, method=rmethod)
		
 		covariatesStack = stack(cellsSmall, covariatesStack)
		covariatesSP = as(covariatesStack, "SpatialPointsDataFrame")
		covariatesDF = covariatesSP@data
		
		notInData = setdiff(names(covariates), names(data))
		covData = stackRasterList(covariates, data, method=rmethod[notInData])
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

	
	callGeneric(
			formula=formula, data=dataDF,
			grid = cellsSmall,
			covariates=covariatesDF, 
			buffer=buffer,
			...
	)
	
	
}

gm.dataSpatial = 
function(formula, data,  grid, 
		covariates=NULL, 
		buffer=0,
		...) {

 

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
			names(data)[apply(data@data,2,is.factor)],
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
		if(!.compareCRS(covariates[[D]], data, unknown=TRUE) ) {
			if(require('rgdal', quietly=TRUE ) ) { 
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


	callGeneric(
			formula=formula, data=data,
			grid = cellsSmall,
			covariates=covariatesDF, 
			...
	)
	
	
}
