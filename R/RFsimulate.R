useRandomFields = function() {
  # false if the RandomFields package shouldnt be used
  requireNamespace("RandomFields", quietly = TRUE) &
    all(options()$useRandomFields)
}

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
  signature(model="RMmodel", x="GridTopology"),
  function(model, x,  data = NULL, 
    err.model=NULL, n = 1, ...)  {
    
    
    if (useRandomFields()) { 
      
      # convert data to an RFspdf (it might be a vanilla spdf)	
      if(!is.null(data)) {
        if(any(class(data)=='SpatialPointsDataFrame')){
          data = RandomFields::conventional2RFspDataFrame(
            data=as.matrix(data@data[,1]), 
            coords=as.matrix(data@coords),
            n = 1)
          # for some reason there is sometimes an error of
          # Error in simu[index, ] : subscript out of bounds
          # in RandomFields unless I do the following
          data@coords[1,1] = data@coords[1,1] + 10e-7
        } }
      
      theArgs = list(...)
      theArgs$x = x
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
      theArgs$spConform=TRUE
      res = do.call(RandomFields::RFsimulate, theArgs)
      res = SpatialGridDataFrame(
        grid = res@grid, data=res@data
      )
      
    } else {
      warning("RandomFields package not available")
      res = NULL
    }
    res
  }
)

setMethod("RFsimulate", 
  signature(model="RMmodel", x="SpatialPoints"),
  function(model, x, 	data = NULL, 
    err.model=NULL, n = 1, ...)  {
    
    if (useRandomFields()) { 
      # convert data to an RFspdf (it might be a vanilla spdf)	
      
      if(!is.null(data)) {
        if(any(class(data)=='SpatialPointsDataFrame')){
          data = RandomFields::conventional2RFspDataFrame(
            array(
              as.matrix(data@data),
              c(dim(data@data),1)
            ),
            coordinates(data),
            n=dim(data)[2]
          )
#		data=as(data, "RFspatialPointsDataFrame") 
        }
      }
      
      theArgs = list(...)
      theArgs$model = model
      theArgs$x = x@coords[,1]
      theArgs$y = x@coords[,2]
      if(!is.null(data))
        theArgs$data = data
      if(!is.null(err.model))
        theArgs$err.model = err.model
      theArgs$n = n
      theArgs$spConform=TRUE
      
      res= try(do.call(RandomFields::RFsimulate, theArgs))
      
    } else { # RandomFields not available
      warning("RandomFields package not available")
      res = NULL
    }
    res
  }
)


setMethod("RFsimulate", 
  signature(model="numeric", x="SpatialPoints"), 
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
        theSim = theSim@data
      }
      names(theSim) = gsub("^variable1(\\.n)?","sim", names(theSim))
      
    } else { #RandomFields not available
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
        theSim = theSim + xcov %*% 
          (Linv %*% data.frame(data)[,1])
      }
      theSim = as.data.frame(as.matrix(theSim))
      names(theSim) = paste("sim", 1:n,sep="")
    } # end no RandomFields	
    
    
    
    if(n==1){
      names(theSim) = gsub("1$", "", names(theSim)) 
    }
    
    res = SpatialPointsDataFrame(
      SpatialPoints(x),
      data=theSim) 
    
    res
  }
)

setMethod("RFsimulate", 
  signature(model="numeric", 
    x="GridTopology"), 
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
        theSim=as.data.frame(matrix(NA, prod(x@cells.dim), n))
      } else {
        theSim = theSim@data
      }
    } else { #RandomFields not available
      modelFull = fillParam(model)
      res2 = raster(x)
      if(ncell(res2) > 500) {
        message("install the RandomFields package for faster simulations")
      }
      if(n>1) {
        res2 = brick(res2, nl=n)
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
        
        
        xRaster = raster(x)
        
        resC = .C(
          C_maternRaster, #"maternRasterConditional",
          # raster
          as.double(xmin(xRaster)), 
            as.double(xres(xRaster)), 
            as.integer(ncol(xRaster)), 
            # raster y
            as.double(ymax(xRaster)), 
            as.double(yres(xRaster)), 
            as.integer(nrow(xRaster)),
            # data 
# as.double(y@data[,1]),    
#          as.double(y@coords[,1]), 
#          as.double(y@coords[,2]), 
#          N=as.integer(Ny), 
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
        
        theCov =  new("dsyMatrix", 
          Dim = as.integer(rep(ncell(xRaster),2)), 
          uplo="L",
          x=resC$result)
        
        Linv = matern(data, param=paramForData, type='inverseCholesky')
        covpreddata = matern(res2, y=data, param=model)
        xcov =  tcrossprod(covpreddata,Linv)
        
        theCov  =theCov - tcrossprod(xcov)
        theChol = chol(theCov)
        
      } else {
        theChol = matern(res2, param=model,
          type='cholesky')
        # do the multiplication with random normals in C?
      }
      
      theSim = crossprod(theChol , theRandom)
      
      if(!is.null(data)) {
        theSim = theSim + xcov %*% 
          (Linv %*% data.frame(data)[,1])	
      }
      theSim = as.data.frame(as.matrix(theSim))
      
      names(theSim) = paste("sim", 1:ncol(theSim),sep="")
    }
    
    if(n==1){
      names(theSim) = gsub("1$", "", names(theSim)) 
    }
    
    
    res = SpatialGridDataFrame(x,theSim)
    
    res
  }
)

setMethod("RFsimulate", 
  signature("matrix", "Raster"),
  function(model, x, data=NULL, 
    err.model=NULL, n = nrow(model), ...)  {
    
    if(is.null(rownames(model)))
      rownames(model) = paste("par", 1:nrow(model),sep="") 
    Siter = round(seq(1,nrow(model), len=n))
    
    if(!is.null(data)) {
      if(!any(class(data)== "SpatialPointsDataFrame"))
        warning("data should be a SpatialPointsDataFrame")
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
      result = brick( 
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
  signature(model="ANY", x="Raster"), 
  function(
    model, x,
    data=NULL, 
    err.model=NULL, n = 1, ...)  {
    
    theproj = x@crs
    x = as(x, "GridTopology")
    res = callGeneric( 
      model=model, x=x,  
      data=data, 
      err.model=err.model, 
      n=n, ... 
    )
    res2 = brick(res)
    if(nlayers(res2)==1){
      res2 = res2[[1]]			
    }
    res2@crs = theproj
    names(res2) = gsub("^variable1\\.n","sim", names(res2))
    
    return(res2)
  }
)


RFsimulate.SPgrid	=	function(
  model, x, 
  data=NULL, 
  err.model=NULL, n = 1, ...)  {
  xOrig = x
  x= as(x, "GridTopology")
  
#			res2 = callGeneric( 
#					model, x,  
#					data=data,	
#					err.model= err.model, 
#          n=n, ... 
#			)
  res2 = callGeneric()
  if(!length(grep("DataFrame$", class(xOrig)))) {
    xOrig = as(xOrig, paste(class(xOrig), "DataFrame",sep=""))
  }
  xOrig@data = res2@data
  
  return(xOrig)
}


setMethod("RFsimulate", 
  signature("numeric", "SpatialPixels"), 
  RFsimulate.SPgrid)


setMethod("RFsimulate", 
  signature("numeric", "SpatialGrid"), 
  RFsimulate.SPgrid)


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
  signature("matrix", "Spatial"),
  function(model, x,  data=NULL, 
    err.model=NULL, n = nrow(model), ...)  {
    
    if(is.null(rownames(model)))
      rownames(model) = paste("par", 1:nrow(model),sep="") 
    
    Siter = round(seq(from=1, to=nrow(model), len=n))
    
    if(!is.null(data)) {
      if(!any(class(data)== "SpatialPointsDataFrame"))
        warning("data should be a SpatialPointsDataFrame")
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
      
      resHere = 	callGeneric(
        model=model[D,], 
        x, data= data[,D], 
        err.model= err.model[D], n=1, ...)
      
      result = cbind(
        resHere@data[,'sim'],
        result
      )
    }
    
    if(n>1)
      colnames(result) = paste('sim', 1:n, sep='')
    
    resHere@data	= as.data.frame(result)
    
    resHere
  }
)


