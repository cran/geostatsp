useRandomFields = function() {
  # false if the RandomFields package shouldnt be used
  requireNamespace("RandomFields", quietly = TRUE) &
    all(options()$useRandomFields)
}

setClass("GridTopology",
  representation(
    cellcentre.offset = 'numeric',
    cellsize = 'numeric',
    cells.dim = 'integer'
  )
)
setClass('RMmodel', 
  representation(
    # call='RMexp(var=1, sclae=1, Aniso=id, proj=id)
    call = "language",
    # name='RMexp'
    name = "character",
    # submodels=NULL, submodels=list(RMmodel1, RMmodel2)
    submodels = "list",
    # model specific parameter 
    par.model = "list",
    # var=1, scale=1, Aniso=id, proj=id 
    par.general = "list"
  )
)


setGeneric('RFsimulate', function(
    model,x, data=NULL, 
    err.model=NULL, n=1, ...) 
    standardGeneric("RFsimulate")
)


setMethod("RFsimulate", 
  signature(model="RMmodel", x="SpatRaster"),
  function(model, x,  data = NULL, 
    err.model=NULL, n = 1, ...)  {
    
    
    if (useRandomFields()) {
      
      # convert data to an RFspdf (it might be a vanilla spdf)	
      if(!is.null(data)) {
        if(any(class(data)=='SpatVector')){
          data = RandomFields::conventional2RFspDataFrame(
            data=as.matrix(values(data)[,1]), 
            coords=crds(data),
            n = 1)
          # for some reason there is sometimes an error of
          # Error in simu[index, ] : subscript out of bounds
          # in RandomFields unless I do the following
#          data@coords[1,1] = data@coords[1,1] + 10e-7
        } }
      
      theArgs = list(...)
      theArgs$x = new("GridTopology",
          cellcentre.offset = c(0,0),
          cellsize = terra::res(x)[1:2],
          cells.dim = as.integer(dim(x)[1:2]))
      if(!is.null(data))
        theArgs$data = data
      if(!is.null(err.model)) {
        if(is.numeric(err.model))
          err.model = RandomFields::RMnugget(
            var=err.model)
        theArgs$err.model = err.model
      }
      
      theArgs$model = model
      theArgs$n = n
      theArgs$spConform=FALSE
      resFromRF = do.call(RandomFields::RFsimulate, theArgs)
      res = rast(x, nlyrs = n)
      values(res) = resFromRF
    } else {
      warning("RandomFields package not available")
      res = rast(x, nlyrs = n)
      values(res) = NA
    }
    res
  }
)

setMethod("RFsimulate", 
  signature(model="RMmodel", x="SpatVector"),
  function(model, x, 	data = NULL, 
    err.model=NULL, n = 1, ...)  {
    
    if (useRandomFields()) { 
      # convert data to an RFspdf (it might be a vanilla spdf)	
      
      if(!is.null(data)) {
        if(any(class(data)=='SpatVector')){
          data = RandomFields::conventional2RFspDataFrame(
            as.matrix(values(data)),
            crds(data),
            n=dim(data)[2]
          )
#		data=as(data, "RFspatialPointsDataFrame") 
        }
      }
      
      xCoords = crds(x)
      theArgs = list(...)
      theArgs$model = model
      theArgs$x = xCoords[,1]
      theArgs$y = xCoords[,2]
      if(!is.null(data))
        theArgs$data = data
      if(!is.null(err.model))
        theArgs$err.model = err.model
      theArgs$n = n
      theArgs$spConform=TRUE
      theArgs$grid = FALSE
      
      res= try(do.call(RandomFields::RFsimulate, theArgs))
      
    } else { # RandomFields not available
      warning("RandomFields package not available")
      res = NULL
    }
    res
  }
)


setMethod("RFsimulate", 
  signature(model="numeric", x="SpatVector"), 
  function(model, x,  data = NULL, 
    err.model=NULL, n = 1, ...)  {
    
    if (useRandomFields()) {
      model = modelRandomFields(model)
      if(!is.null(err.model))
        err.model = RandomFields::RMnugget(var=err.model)
      
      theSim=callGeneric(
        model, x,  data = err.model, 
        err.model=err.model, n = n, ...
      )
      #model, x,  
#		data=data,	err.model= err.model, n=n  ,  ...)
      if(all(class(theSim)%in%c('try-error', 'NULL'))) {
        warning("error in RandomFields")
        theSim=as.data.frame(matrix(NA, length(x), n))
      } else {
        theSim = values(theSim)
      }
      names(theSim) = gsub("^variable1(\\.n)?","sim", names(theSim))
      
    } else { #RandomFields not available
      model = fillParam(model)
      model['nugget']=0
      if(!is.null(data)) {
        theCov = matern(x, param=model)
        #covd = matern(data, param=model)
        #if(!is.null(err.model))
        #	diag(covd) = diag(covd) + err.model
        #Linv = solve(chol(covd))
        paramForData = model
        if(!is.null(err.model))
          paramForData['nugget']=as.numeric(err.model)
        Linv = matern(data, param=paramForData,
          type='inverseCholesky')
        covpreddata = matern(x, y=data, param=model)
        xcov =  tcrossprod( covpreddata,Linv)
        theCov  =theCov - tcrossprod(xcov)
        theChol = chol(theCov)
      } else {
        theChol = matern(x, param=model, type='cholesky')
      }
      theRandom = matrix(
        rnorm(n*nrow(theChol)), 
        nrow=nrow(theChol), ncol=n)
      theSim = theChol %*% theRandom
      if(!is.null(data)) {
        LinvData = Linv %*% data.frame(data)[,1]
        LinvDataCov0  = xcov %*% LinvData
        theSim = theSim + LinvDataCov0
      }
      theSim = as.data.frame(as.matrix(theSim))
      names(theSim) = paste("sim", 1:n,sep="")
    } # end no RandomFields	
    
    
    
    if(n==1){
      names(theSim) = gsub("1$", "", names(theSim)) 
    }
    
    res = deepcopy(x)
    terra::values(res) = theSim
    res
  }
)

setMethod("RFsimulate", 
  signature(model="numeric", 
    x="SpatRaster"), 
  function(model, x, data = NULL, 
    err.model=NULL, n = 1, ...)  {
    if (useRandomFields()) { 
      model = modelRandomFields(model)
      if(!is.null(err.model))
        err.model = RandomFields::RMnugget(var=err.model)
      
      theSim=callGeneric(model, x,  
        data=data,	err.model= err.model,
        n=n  , 
        ...)
      
      names(theSim) = gsub("^variable1(\\.n)?","sim", names(theSim))
      
      if(all(class(theSim)%in%c('try-error', 'NULL'))) {
        warning("error in RandomFields")
        theSim=as.data.frame(matrix(NA, prod(dim(x)[1:2], n)))
      } else {
        theSim = values(theSim)
      }
    } else { #RandomFields not available
      modelFull = fillParam(model)
      res2 = rast(x)
      if(ncell(res2) > 500) {
        message("install the RandomFields package for faster simulations")
      }
      if(n>1) {
        res2 = rast(res2, nlyrs=n)
      }
      
      theRandom = matrix(
        rnorm(n*ncell(res2)), 
        nrow=ncell(res2), 
        ncol=n)
      
      if(!is.null(data)) {
        #covd = matern(data, param=model)
        #if(!is.null(err.model))
        #	diag(covd) = diag(covd) + err.model
        #Linv = solve(chol(covd))
#        paramForData = model
#        if(!is.null(err.model))
#          theCov = matern(res2, param=model)
        
        paramForData = fillParam(model)
        if(!is.null(err.model))
          paramForData['nugget']=as.numeric(err.model)
        
        
        #covd = matern(data, param=model)
        
        #if(!is.null(err.model))
        # diag(covd) = diag(covd) + err.model
        #Linv = solve(chol(covd))
        
        # do all this in C?
        
        
        xRaster = rast(x)
        
        resC = .C(
          C_maternRaster, #"maternRasterConditional",
          # raster
          as.double(xmin(xRaster)), 
            as.double(res(xRaster)[1]), 
            as.integer(ncol(xRaster)), 
            # raster y
            as.double(ymax(xRaster)), 
            as.double(res(xRaster)[2]), 
            as.integer(nrow(xRaster)),

            # conditional covariance matrix
            result=as.double(matrix(0, ncell(xRaster), ncell(xRaster))),
            # parameters
            as.double(modelFull["range"]),
            as.double(modelFull["shape"]),
            as.double(modelFull["variance"]),
            as.double(modelFull["anisoRatio"]),
            as.double(modelFull["anisoAngleRadians"]),
            as.integer(1)
          )
        
        theCov =  new("dpoMatrix", 
          Dim = as.integer(rep(ncell(xRaster),2)), 
          uplo="L",
          x=resC$result)
        
        Linv = matern(data, param=paramForData, type='inverseCholesky')
        covpreddata = matern(res2, y=data, param=model)
        xcov =  tcrossprod(covpreddata,Linv)
        xcov2 = tcrossprod(xcov)

        theCov  =theCov - xcov2
        theChol = chol(theCov)
        
      } else { # not conditional simulation
        theChol = matern(res2, param=model,
          type='cholesky')
        # do the multiplication with random normals in C?
      }
      
      theSim = crossprod(theChol , theRandom)
      
      if(!is.null(data)) {
        # bug here?
        LinvData = Linv %*% data.frame(data)[,1]
        LinvDataCov1 = xcov %*% LinvData
        theSim = theSim + matrix(LinvDataCov1, nrow(theSim), ncol(theSim))

      }
      theSim = as.data.frame(as.matrix(theSim))
      
      names(theSim) = paste("sim", 1:ncol(theSim),sep="")
    }
    
    if(n==1){
      names(theSim) = gsub("1$", "", names(theSim)) 
    }
    
    
    result = rast(x, nlyrs = ncol(theSim), vals = theSim, 
      names = colnames(theSim))
  }
)

setMethod("RFsimulate", 
  signature("matrix", "SpatRaster"),
  function(model, x, data=NULL, 
    err.model=NULL, n = nrow(model), ...)  {
    
    if(is.null(rownames(model)))
      rownames(model) = paste("par", 1:nrow(model),sep="") 
    Siter = round(seq(1,nrow(model), len=n))
    
    if(!is.null(data)) {
      if(!any(class(data)== "SpatVector"))
        warning("data should be a SpatVector")
      # check data variables
      if(ncol(data) == 1) {
        data = data[,rep(1,max(Siter))]
        
      } else if(ncol(data) > 1) {
        # if there's more than one, assume we're interested in the first one
        
        if(ncol(data)!= nrow(model)){
          warning("number of columnns in data should be either 1 or equal to number of rows of model")
        }
      }
    } else {
      # do something so data[,D] doesn't break
      data = NULL
    }
    if(!is.null(err.model)) {
      if(is.numeric(err.model)){
        if(length(err.model)==1) {
          err.model = rep(err.model, length(Siter))
        } else if(length(err.model)!=nrow(model)){
          warning("number of values in err.model should be either 1 or equal to number of rows of model")
        }
      }
    } else {
      err.model= NULL
    }
    
    
    result = NULL
    
    for(D in rev(Siter)) {
      resultHere = 	callGeneric(
        model=model[D,], 
        x=x, 
        data= data[,D], 
        err.model=err.model[D], 
        n=1,
        ...)
      result = c( 
        resultHere,result
      )
    }
    
    if(n>1) {
      names(result) = paste('sim', 1:n, sep='')
    } else {
      names(result) = 'sim'
    }
    result
  }
)





setMethod("RFsimulate", 
  signature("data.frame", "ANY"),
  function(model, x, data=NULL, 
    err.model=NULL, n = nrow(model), ...)  {
    
    model = as(model, "matrix")
    
    callGeneric() 
#				model, x, 
#				data=data,	err.model= err.model,
#				n=n  , ... 
#		)
  }
)

setMethod("RFsimulate", 
  signature("matrix", "SpatVector"),
  function(model, x,  data=NULL, 
    err.model=NULL, n = nrow(model), ...)  {
    
    if(is.null(rownames(model)))
      rownames(model) = paste("par", 1:nrow(model),sep="") 
    
    Siter = round(seq(from=1, to=nrow(model), len=n))
    
    if(!is.null(data)) {
      if(!any(class(data)== "SpatVector"))
        warning("data should be a SpatVector")
      # check data variables
      if(ncol(data) == 1) {
        data = data[,rep(1,max(Siter))]
        
      } else if(ncol(data) > 1) {
        # if there's more than one, assume we're interested in the first one
        
        if(ncol(data)!= nrow(model)){
          warning("number of columnns in data should be either 1 or equal to number of rows of model")
        }
        
      }
    } else {
      # do something so data[,D] doesn't break
      data = NULL
    }
    if(!is.null(err.model)) {
      if(is.numeric(err.model)){
        if(length(err.model)==1) {
          err.model = rep(err.model, length(Siter))
        } else if(length(err.model)==nrow(model)){
          err.model = err.model[Siter]
        } else {
          warning("number of values in err.model should be either 1 or equal to number of rows of model")
        }
      }
    } else {
      err.model= NULL
    }
    
    result = NULL
    
    for(D in rev(Siter)) {
      
      resHere = callGeneric(
        model=model[D,], 
        x, data= data[,D], 
        err.model= err.model[D], n=1, ...)
      
      result = cbind(
        values(resHere)[,'sim'],
        result
      )
    }
    
    if(n>1)
      colnames(result) = paste('sim', 1:n, sep='')
    
    values(resHere)	= as.data.frame(result)
    
    resHere
  }
)


