setGeneric('glgm', 
  function(
    formula, 
    data,
    grid,
    covariates,
    buffer=0,  
    shape=1, 
    prior, 
    ...) {
    standardGeneric("glgm")
  }
  )


setMethod("glgm", 
  signature(formula="ANY", data="ANY", grid="ANY", covariates="ANY"), 
  function( 
    formula, 
    data,
    grid,
    covariates,
    buffer=0,  
    shape=1, 
    prior, 
    ...) {

# clean names of covariates
    if(!missing(covariates)) {
      if(is.list(covariates)) {
        if(!length(names(covariates))) {
          names(covariates) = paste0("covariate", 1:length(covariates))
        }
      }
      if(is.list(covariates)|is.data.frame(covariates)) {
        names(covariates) = make.names(names(covariates))
      }
    }


# numeric or missing formula
    if(missing(formula)) formula = 1
    if(is.null(formula)) formula = 1
    if(is.numeric(formula)) {
      if(missing(data)) 
        stop("data must be provided if formula is numeric")
      if(!length(names(data)))
        stop("data must have names if formula is numeric")
      formula = names(data)[formula]
    }

# character formula
    if(is.character(formula)) {
      if(length(formula)==1)
        formula = unique(c(formula, names(covariates)))
      if(length(formula)==1)
        formula = c(formula, '1')

      formula = paste(formula[1] , "~",
        paste(formula[-1], collapse=" + ")
        )
      formula = as.formula(formula)
    }

    callGeneric()
  }
  )







# extrat covariates for data, convert covariates to a stack
setMethod("glgm", 
  signature(formula="formula", data="Raster", grid="ANY", covariates="ANY"),
  function(
    formula, 
    data,
    grid,
    covariates,
    buffer=0,  
    shape=1, 
    prior, 
    ...) {

    if(is.numeric(grid))
      grid = squareRaster(data, grid)

    dataCov = gm.dataRaster(
      formula, data,
      grid,
      covariates,
      buffer)


 #   data = dataCov$data 
 #   grid = dataCov$grid
 #   covariates = dataCov$covariates

   callGeneric(
      formula = dataCov$formula,
      data = dataCov$data,
      grid = dataCov$grid,
      covariates = dataCov$covariates,
      buffer,
      shape,
      prior, 
      ...
      )
   }
  )


setMethod("glgm", 
  signature(formula="formula", data="Spatial", grid="ANY", covariates="ANY"),
  function(    formula, 
    data,
    grid,
    covariates,
    buffer=0,  
    shape=1, 
    prior, 
    ...) {

    if(is.numeric(grid))
      grid = squareRaster(data, grid)

    dataCov = gm.dataSpatial(
      formula, data, 
      grid, covariates, buffer)

#    data = dataCov$data@data 
#    grid = dataCov$grid
#    covariates = dataCov$covariates

    callGeneric(
      formula,
      dataCov$data@data,
      dataCov$grid,
      dataCov$covariates,
      buffer,
      shape,
      prior, 
      ...
      )
  }
  )

#################
#### the real work
##################

setMethod("glgm", 
  signature(
    formula="formula", 
    data="data.frame",
    grid="Raster", 
    covariates="data.frame"), 
  function(
    formula, 
    data,
    grid,
    covariates,
    buffer=0,  
    shape=1, 
    prior, 
    ...) {

# undocumented options for ...
    getRidDots = c(
      'priorCI', # legacy prior specification
      'spaceExtra', # extra arguments for the space formula
      'spaceFormula') # override the space formula

    if(!any(names(grid)=='space')) {
      grid = setValues(raster(grid), 1:ncell(grid))
      names(grid) = 'space'
    }

    allVars = allVarsP(formula)

    if(!all(allVars %in% names(data)) )
      warning("some covariates seem to be missing: formula ", 
        paste(allVars, collapse=" "), ", data: ", 
        paste(names(data), collapse=" "))

    cells = raster::trim(grid[['space']])
    firstCell = values(cells)[1]
    cellDim = dim(cells)[1:2]
    # first cell = 2 * buffer^2 + ncolSmall * buffer + buffer
    # buffer = -(nrowSmall+1) + sqrt (  (nrowSmall+1)^2 + 8 firstCell / 4
    buffer = (-(cellDim[1]+1) + sqrt(  (cellDim[1]+1)^2 + 8* (firstCell-1) ))/4
    # data, cells, and covariates must have varilable called 'space'		
    # values of cells must be index numbers, and cells shouldnt include the buffer		
    forInla = thedots = list(...)
    forInla = forInla[setdiff(names(forInla), getRidDots)]

    if(!any(names(thedots)=="family")) {
      forInla$family =  "gaussian"
    }

    if(any(names(thedots)== 'priorCI')) {
    # legacy priors
      priorList = priorLegacy(thedots$priorCI, forInla$family, cellSize = xres(cells))
    } else {
      if(missing(prior)) prior=list(range=NULL)
        priorList = priorInla(prior, forInla$family, cellSize = xres(cells))
    }

    # done priors

    # formula for spatial random effect
    spaceFormula = c(
      ".~.+ f(space, model='matern2d', ",
      "nrow=", nrow(cells)+2*buffer, 
      ", ncol=", ncol(cells)+2*buffer,
      ", nu=", shape, 
      ", hyper = list(",
      "range=",
      priorList$range$string,
      ",",
      "prec=",
      priorList$sd$string,
      " )", 
      paste(c('', thedots$spaceExtra), collapse=', '),
      " )", 
      sep=""
      )

    formulaOrig = formula

    if(any(names(thedots)=='spaceFormula')) {
      formula = update.formula(formulaOrig,	as.formula(thedots$spaceFormula))
    } else {
      formula = update.formula(formulaOrig,	as.formula(paste(spaceFormula, collapse=' ')))
    }


    # sort out factors
    # variables with 'factor(x)' in the formula
    thevars = c(allVarsP(formulaOrig), 'space')
    thevars = grep("^factor\\(", thevars, value=TRUE)

    varsInData = unlist(lapply(data, is.factor))
    varsInData = names(data)[varsInData]
    varsInData = intersect(all.vars(formula), varsInData)

    thevars = union(varsInData, thevars)

    if(length(thevars)){

      thevars = gsub("^factor\\(|\\)", "", thevars)
      # loop through factors
      for(D in thevars){
        # biggest category is baseline
        thetable = table(data[,D])
        thebase = names(sort(thetable,decreasing=TRUE))[1]
        newLevels = unique(c(thebase, levels(factor(data[,D]))))
        data[,D] = factor(data[,D], levels=newLevels)
        if(D %in% colnames(covariates))
          covariates[,D] = factor(
            covariates[,D],
            levels=levels(data[,D]))
      }
    }

    # variables wrapped in an f()
    termsInF = grep("^f[(]", attributes(allVarsP(formulaOrig))$orig, value=TRUE)
    termsInF = gsub("^f[(]|[,].*|[[:space:]]", "", termsInF)
    termsInF = intersect(termsInF, colnames(covariates))

    theFactors = thevars

    if(!any(names(thedots)=='lincomb')) {
      # create linear combinations object for prediction.
      # create formula, strip out left variable and f(...) terms
      formulaForLincombs = unlist(strsplit(as.character(formulaOrig), "~"))
      formulaForLincombs = formulaForLincombs[length(formulaForLincombs)]

#          rownames(attributes(terms(formulaOrig))$factors), 
#          value=TRUE)
      covariatesInF = covariates[,termsInF, drop=FALSE]

      formulaForLincombs =
      gsub("\\+?[[:space:]]*f\\(.*\\)[[:space:]]?($|\\+)", "+", 
        formulaForLincombs)
      # strip out offsets
      formulaForLincombs =
      gsub("\\+?[[:space:]]*offset\\([[:print:]]*\\)[[:space:]]?($|\\+)", "+", formulaForLincombs)

      # convert multiple + to a single +
      formulaForLincombs = gsub(
        "\\+[[:space:]]?\\+([[:space:]]?\\+)?", "+",
        formulaForLincombs)
      # strip out trailing +
      formulaForLincombs = gsub("\\+[[:space:]]?$", "", formulaForLincombs)

      # if we have covariates in the formula and in the data
      if(max(nchar(formulaForLincombs)) & nrow(covariates) ) {

        formulaForLincombs=as.formula(paste("~", formulaForLincombs))

        # variables in the model but not in prediction rasters
        thevars = allVarsP(formulaForLincombs) #rownames(attributes(terms(formulaForLincombs))$factors)
        thevars = gsub("^factor\\(|\\)", "", thevars)
        varsInPredict = thevars[thevars %in% names(covariates)]
        cantPredict = thevars[! thevars %in% names(covariates)]
        theFactors2 = grep("factor\\([[:print:]]*\\)", cantPredict)
        if(length(theFactors2)) {
          temp = cantPredict
          cantPredict = cantPredict[-theFactors2]
          theFactorsInFormula = temp[theFactors2]
        }
        if(length(cantPredict)){
          covariates[,cantPredict]= 0
        }
        covariates = covariates[,c("space", thevars),drop=FALSE]
        lincombMat = model.matrix(update.formula(
          formulaForLincombs, ~.+space),
        covariates, na.action=NULL)
        lincombMat[lincombMat==0] = NA

        spaceCol = grep("^space$", colnames(lincombMat), value=TRUE, ignore.case=TRUE)

        thelincombs <- apply(lincombMat, 1, lcOneRow, idxCol=spaceCol)
        names(thelincombs) = paste("c", lincombMat[,spaceCol],sep="")

      } else { # no covariates or no INLA
      thelincombs=list()	
      for(D in 1:ncell(cells)) {
        thelincombs[[D]] = list(
          list(
            "(Intercept)"=list(weight=1)
            ),
          list(
            space=list(weight=1, idx=values(cells)[D])
            )
          )
      }
      names(thelincombs) = paste("c", values(cells),sep="")
    }



      # add in the covariates wrapped in f
    for(Dcov in termsInF) {
      for(Dtext in names(thelincombs)) {
        Dnum = as.numeric(gsub("c", '', Dtext))
        Dcell = match(Dnum, covariates[,'space'])
        if(!is.na(covariatesInF[Dcell,Dcov])) {
          covToAdd = list(list(
            weight = 1,
            idx = covariatesInF[Dcell,Dcov]
            ))
          names(covToAdd) = Dcov

          thelincombs[[Dtext]] = c(
            thelincombs[[Dtext]],
            list(
              covToAdd
              )
            ) 
        } else {
          thelincombs[[Dtext]] = NULL
        }
      }
    }

    getRid = which(unlist(lapply(thelincombs, is.null)))
    if(length(getRid)) thelincombs = thelincombs[-getRid]

    forInla$lincomb = thelincombs
  } # end adding lincombs unless lincomb is user-supplied


    # get rid of observations with NA's in covariates
  allVars = allVarsP(formulaOrig)

  if(length(allVars)) {
    theNA = apply(data[,c(allVars,'space'),drop=FALSE], 
      1, function(qq) any(is.na(qq)))
  } else {
    theNA = rep(FALSE, ncol(data))
  }

  data = data[!theNA,]
  forInla$data = data
  forInla$formula = formula
  if(any(names(thedots)=='Ntrials'))
    forInla$Ntrials = thedots$Ntrials[!theNA]
  if(any(names(thedots)=='weights'))
    forInla$weights = thedots$weights[!theNA]



    # if model is gaussian, add prior for nugget
  if(!is.null(priorList$sdObs)) {
    if(!length(forInla$control.family$hyper$prec)) {
      forInla$control.family$hyper$prec =
      eval(parse(text=priorList$sdObs$string))
    } else {
      priorList = priorList[setdiff(names(priorList), 'sdObs')]
    }
  }

  familyShapeName = grep("familyShape", names(priorList), value=TRUE)
  if(length(familyShapeName)) {
    forInla$control.family$hyper$theta =
    priorList$familyShapePrior 
  }


    # get rid of some elements of forInla that aren't required
  if(!length(forInla$lincomb)) 
    forInla = forInla[setdiff(names(forInla), 'lincomb')] 


  if(requireNamespace("INLA", quietly=TRUE)) {
    if(identical(forInla$verbose, TRUE)) {
      tFile = tempfile('glgm', tempdir(), '.rds')
      message(paste('saving INLA objects as', tFile))
      saveRDS(forInla, file=tFile)
    }
    inlaResult = try(do.call(INLA::inla, forInla))
  } else {
    inlaResult = 
    list(logfile="INLA is not installed. \n see www.r-inla.org")
  }
  if(identical(forInla$verbose, TRUE)) {
    message("inla done") 
  }
  if(all(names(inlaResult)=="logfile") | class(inlaResult) == 'try-error')
    return(c(forInla, list(inlares=inlaResult, prior = priorList)))


  params = list(range=list(), scale=list())

    # posterior distributions for range

  # convert from cells to spatial units (i.e. km)
  if("Range for space" %in% names(inlaResult$marginals.hyperpar)) {
  params$range$posterior = 
    inlaResult$marginals.hyperpar[["Range for space"]] %*% 
    diag(c(xres(cells), 1/xres(cells)))
  params$range$posterior = cbind(
    x=params$range$posterior[,1],
    y=params$range$posterior[,2],
    prior=priorList$range$dprior$range(params$range$posterior[,1])
    )

  params$range$prior = priorList$range[setdiff(names(priorList$range), 'dprior')]
  params$range$dprior = priorList$range$dprior$range
  params$scale$dprior = priorList$range$dprior$scale

  params$scale$posterior = cbind(
    x = 1/params$range$posterior[,'x'],
    y = params$range$posterior[,'y'] * params$range$posterior[,'x']^(2)
    )


  theMaxX = max(params$scale$posterior[
   params$scale$posterior[,'y'] > 10^(-4) * max(params$scale$posterior[,'y']), 'x'])

  params$scale$posterior = stats::approx(
    params$scale$posterior[,'x'],   
    params$scale$posterior[,'y'], 
    seq(0, theMaxX, len=nrow(params$scale$posterior)),
    rule = 2
    )

  params$scale$posterior = cbind(
    do.call(cbind, params$scale$posterior),
    prior = priorList$range$dprior$scale(params$scale$posterior$x)
    )

# get rid of right tail
  params$range$posterior = params$range$posterior[
  params$range$posterior[,'y'] > 10^(-4) * max(params$range$posterior[,'y']) |
  params$range$posterior[,'x'] < 10^(-1) * max(params$range$posterior[,'x']), ]


  params$range$postK = params$range$posterior %*% diag(c(1/1000, 1000, 1000))
  params$scale$postK = params$scale$posterior %*% diag(c(1000, 1/1000, 1/1000))
  colnames(params$range$postK) = colnames(params$scale$postK) = 
  colnames(params$range$posterior)
}
  Ssd =grep("^sd", names(priorList), value=TRUE)
  names(Ssd)[grep("^sd$", names(priorList))] = 'Precision for space'
  names(Ssd)[grep("^sdObs$", names(priorList))] = 'Precision for the Gaussian observations'

  for(Dsd in names(Ssd)) {
    DsdName = Ssd[Dsd]
    params[[DsdName]]= list(
      posterior = precToSd(inlaResult$marginals.hyperpar[[Dsd]]),
      prior = priorList[[DsdName]])
    if(!is.null(priorList[[DsdName]]$dprior))
      params[[DsdName]]$posterior = cbind(
        params[[DsdName]]$posterior,
        prior = priorList[[DsdName]]$dprior(params[[DsdName]]$posterior[,'x']))
  }

  if(length(familyShapeName)) {
    params[[familyShapeName]] = list(
      priorList[[familyShapeName]]$prior
      )
  }

    # random into raster
# E exp(random)

  if("summary.random" %in% names(inlaResult)) {

    temp=unlist(
      lapply(inlaResult$marginals.random$space, function(qq) {
        sum(
          exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
          )
      })
      )
    inlaResult$summary.random[['space']][,"exp"] = temp



    forRast = 	as.matrix(inlaResult$summary.random[["space"]][values(cells),])
    resRasterRandom = 
    brick(extent(cells), nrows=nrow(cells),
      ncols=ncol(cells), crs=projection(cells),
      nl=dim(forRast)[2])
    names(resRasterRandom) = 
    paste("random.", colnames(forRast),sep="")

    values(resRasterRandom) = as.vector(forRast)

  } else {
    return(list(inla=inlaResult, parameters=params))
  }

  inlaResult$marginals.random$space = inlaResult$marginals.random$space[values(cells)]

    # put linear combinations into the raster
    # E exp(lincombs)
  temp=unlist(
    lapply(inlaResult$marginals.lincomb.derived, function(qq) {
      sum(
        exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
        )
    })
    )
  inlaResult$summary.lincomb.derived[,"exp"] = temp

    # E inv logit(lincombs)
  if(length(grep("logit",inlaResult$misc$linkfunctions$names))) {
    temp=unlist(
      lapply(inlaResult$marginals.lincomb.derived, function(qq) {
        eqqx = exp(qq[,"x"])
        sum(
          eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
          )
      })
      )
    inlaResult$summary.lincomb.derived[,"invlogit"] = temp		
  }


    # lincombs into raster
    # theSpaceName will be empty if lincomb is user-supplied
  theSpaceName = grep("^c[[:digit:]]+$", names(inlaResult$marginals.lincomb.derived), value=TRUE)
  theSpace = as.integer(gsub("^c", "", theSpaceName))

  if(length(theSpaceName)) {
    linc = inlaResult$summary.lincomb.derived[theSpaceName,]
    linc$space = theSpace
    inlaResult$marginals.predict = 
    inlaResult$marginals.lincomb.derived

    missingCells = values(cells)[! values(cells) %in% theSpace]

    if(length(missingCells)) {
      toadd = matrix(NA, length(missingCells), dim(linc)[2], 
        dimnames=list(
          paste("c", missingCells, sep=""), 
          colnames(linc)
          )
        )
      toadd[,"space"] = missingCells

      linc = rbind(linc, toadd)

        # Add in empty lists for the marginals of missing cells

      missingMarginals = vector("list", length(missingCells))
      names(missingMarginals) = rownames(toadd)

      inlaResult$marginals.predict = c(	
        inlaResult$marginals.predict,
        missingMarginals)

    }
    linc = as.matrix(linc[match( values(cells), linc$space),])
    inlaResult$marginals.predict = 
    inlaResult$marginals.predict[
    paste("c", values(cells), sep="")
    ]


    resRasterFitted = 
    brick(extent(cells), nrows=nrow(cells),
      ncols=ncol(cells), crs=projection(cells),
      nl=ncol(linc))
    names(resRasterFitted) = 
    paste("predict.", colnames(linc),sep="")

    values(resRasterFitted) = as.vector(linc)
  }  else {
    resRasterFitted = NULL
  }    




# sum(c(0,diff(params$range$posterior[,"x"])) * params$range$posterior[,"y"])
# sum(c(0,diff(params$range$prior[,"x"])) * params$range$prior[,"y"])


  params$summary = inlaResult$summary.fixed

  params$summary = cbind(params$summary, 
    meanExp = unlist(
      lapply(inlaResult$marginals.fixed,
        function(qq) {
          sum(
            exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
            )
        }
        ))
    )

  if(length(grep("logit",inlaResult$misc$linkfunctions$names))) {
    params$summary = cbind(params$summary, 
      meanInvLogit = unlist(
        lapply(inlaResult$marginals.fixed, function(qq) {
          eqqx = exp(qq[,"x"])
          sum(
            eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
            )
        }
        )
        ))
  }


  thecols = paste(c("0.975", "0.5","0.025"), "quant", sep="")

  thesd = c(
    sdNugget= grep("^Precision[[:print:]]*Gaussian observations$", 
      names(inlaResult$marginals.hyperpar), value=TRUE),
    gammaShape = grep("^Precision[[:print:]]*Gamma observations$", 
      names(inlaResult$marginals.hyperpar), value=TRUE),
    weibullShape = grep("^alpha[[:print:]]*weibull", 
      names(inlaResult$marginals.hyperpar), value=TRUE),
    sd = grep("^Precision[[:print:]]*space$", 
      names(inlaResult$marginals.hyperpar), value=TRUE)
    )

  toAddSd = setdiff(
    grep('^Precision', names(inlaResult$marginals.hyperpar), 
      value = TRUE), thesd)
  if(length(toAddSd)) {
    names(toAddSd) = gsub("^Precision for ", "sd ", toAddSd)
    thesd = c(
      thesd,   
      toAddSd
      )
  }      


  params$summary = rbind(params$summary,
    matrix(NA, nrow=length(thesd)+1, ncol=ncol(params$summary),
      dimnames = list(c("range", names(thesd)), 
        colnames(params$summary)))
    )

  names(inlaResult$all.hyper[['random']]) = 
  unlist(lapply(inlaResult$all.hyper[['random']], 
    function(qq) qq$hyperid))

# convert precisions to standard deviations
  for(Dsd in grep("Shape$", names(thesd), invert=TRUE, value=TRUE)) {

    params$summary[Dsd, thecols] = 
    1/sqrt(inlaResult$summary.hyperpar[
      thesd[Dsd],rev(thecols)])
    params$summary[Dsd,"mode"] = 
    1/sqrt(inlaResult$summary.hyperpar[
      thesd[Dsd],'mode'])


    params$summary[Dsd,"mean"] =sum(
      1/sqrt(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"])*
      c(0,diff(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"]))*
      inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"y"]
      )
    params$summary[Dsd,"sd"] =sum(
      1/inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"]*
      c(0,diff(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"]))*
      inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"y"]
      ) - params$summary[Dsd,"mean"]^2

  }

    # shape parameters (gamma or weibull)
  for(Dsd in grep("Shape$", names(thesd), value=TRUE) ) {
    params[[Dsd]]$posterior=
    inlaResult$marginals.hyperpar[[thesd[Dsd]]]

    params$summary[Dsd, colnames(inlaResult$summary.hyperpar)] = 
    inlaResult$summary.hyperpar[thesd[Dsd],]

      # prior
      # if Dsd matches the family argument
    if(inlaResult$all.hyper$family[[1]]$label == gsub("Shape$", '', Dsd)) {
      fPrior = inlaResult$all.hyper$family[[1]]$hyper$theta
      fDist = fPrior$prior
      fParam = fPrior$param
      if(fDist == 'loggamma') {
        xSeq = c(
          stats::qgamma(c(0.001, 0.999), shape=fParam[1], rate=fParam[2]),
          range(params[[Dsd]]$posterior[,1]))
        xSeq = sort(unique(c(0,
          seq(min(xSeq), max(xSeq), len=999)))
        )
        params[[Dsd]]$prior = cbind(
          x = xSeq,
          y = stats::dgamma(xSeq, shape=fParam[1], rate=fParam[2])
          )
      } else if (fDist == 'normal') {
        xSeq = c(
          stats::qlnorm(c(0.1, 0.9), meanlog=fParam[1], sdlog=1/sqrt(fParam[2])),
          range(params[[Dsd]]$posterior[,1]))
        xSeq = sort(unique(c(0,
          seq(min(xSeq), max(xSeq), len=999)))
        )
        params[[Dsd]]$prior = cbind(
          x = xSeq,
          y = stats::dlnorm(xSeq, meanlog=fParam[1], sdlog=1/sqrt(fParam[2]))
          )
      }

    }

  }

  # any additional standard deviations
  for(Dsd in names(toAddSd)) {
    Dprec = toAddSd[Dsd]
    DsdVar = gsub("^sd[[:space:]]+", "", Dsd)

    params[[Dsd]]= list(
      posterior = precToSd(inlaResult$marginals.hyperpar[[Dprec]]),
      prior = 
        inlaResult$all.hyper$random[[DsdVar]]$hyper$theta[c('prior','param')]
        )

    if(!is.null(params[[Dsd]]$prior$prior)) {
      if(params[[Dsd]]$prior$prior =='pc.prec') {
        params[[Dsd]]$prior = priorInla(list(sd = params[[Dsd]]$prior$param))$sd

        params[[Dsd]]$posterior = cbind(
          params[[Dsd]]$posterior,
          prior = params[[Dsd]]$prior$dprior(
            params[[Dsd]]$posterior[,'x'])
          )
      }
    } # end not null prior
  }




# put range in summary, in units of distance, not numbers of cells
  thecolsFull =c("mean","sd",thecols,"mode") 
  params$summary["range",thecolsFull]=				
  xres(cells)*
  inlaResult$summary.hyperpar[
  "Range for space",
  thecolsFull
  ]
  if(any(params$summary["range", 'mean'] > 1000, na.rm=TRUE)) {

    params$summary["range",thecolsFull]=				
    params$summary["range",thecolsFull] / 1000
    rownames(params$summary) = gsub(
      "^range$", "range/1000",
      rownames(params$summary)
      )
  }
  dimnames(params$summary) = lapply(dimnames(params$summary),
    function(qq) {
      qq=gsub("_", " ", qq)
      qq=gsub("\\$", " ", qq)
      qq=gsub("<", " lt ", qq)
      qq=gsub(">", " gt ", qq)
      qq
    }
    )
  params$summary = as.data.frame(params$summary)

  for(Dvar in names(covariates)) {
    theLevels =levels(covariates[[Dvar]])[[1]]
    if(!is.null(nrow(theLevels))){
      for(D in 1:nrow(theLevels)) {
        rownames(params$summary) = gsub(
          paste("(factor)?(\\()?", Dvar, "(\\))?:?", 
            theLevels[D,1],"$",sep=""),
          paste(Dvar, ":",theLevels[D,2],sep=""), 
          rownames(params$summary))
      }
    }
  }

getRid = c('random.ID', 'predict.ID', 'predict.space', 'predict.kld')

  resRaster=brick(stack(
    resRasterRandom[[setdiff(names(resRasterRandom), getRid)]], 
    resRasterFitted[[setdiff(names(resRasterFitted), getRid)]], 
    cells[[setdiff(names(cells), getRid)]]))

  result=list(inla=inlaResult,
    raster=resRaster,
    parameters=params
    )

  result

}
)
