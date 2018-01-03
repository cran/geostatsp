lcOneRow = function(thisrow, idxCol=NULL) {
  thisrow = thisrow[!is.na(thisrow)]
  if(length(thisrow)) {
    thisrow = sapply(thisrow, function(qq) list(list(weight=qq)))
    for(D  in idxCol)
      thisrow[[D]] = list(
        weight=1, 
        idx=thisrow[[D]]$weight
        )
    for(D in names(thisrow))
      thisrow[[D]] = thisrow[D]
    names(thisrow) = paste("v", 1:length(thisrow), sep="")
  }
  thisrow
}

setGeneric('glgm', 
  function(
    formula, data, grid, 
    covariates=NULL, 
    ...) {
    standardGeneric("glgm")
  }
  )

# sort out formula
# null formula
setMethod("glgm", 
  signature("NULL"), 
  function(formula, data, grid, 
    covariates=NULL, ...) {
    formula =  1 
    callGeneric(formula, data, grid, covariates, ...)
  }
  )


setMethod("glgm", 
  signature("numeric"),  
  function(formula, data, grid, 
    covariates=NULL, ...) {

    formula = names(data)[formula]
    callGeneric(formula, data, grid, covariates, ...)
  }
  )

# change character to formula
setMethod("glgm", 
  signature("character"),  
  function(formula, data, grid, 
    covariates=NULL, ...) {

    if(length(names(covariates)))
      names(covariates) = gsub("[[:punct:]]|[[:space:]]","_", names(covariates))
    if(length(covariates) & !length(names(covariates))) 
      names(covariates) = paste("c", 1:length(covariates),sep="")			

    if(length(formula)==1)
      formula = unique(c(formula, names(covariates)))
    if(length(formula)==1)
      formula = c(formula, '1')

    formula = paste(formula[1] , "~",
      paste(formula[-1], collapse=" + ")
      )
    formula = as.formula(formula)

    callGeneric(formula, data, grid, covariates, ...)
  }
  )


# numeric cells, create raster from data bounding box

setMethod("glgm", 
  signature("formula", "ANY", "numeric", "ANY"),
  function(formula, data, grid, covariates=NULL, ...) {
    grid = squareRaster(data, grid)
    callGeneric(formula, data, grid, covariates, ...)
  }
  )



# extrat covariates for data, convert covariates to a stack
setMethod("glgm", 
  signature("formula", "Raster", "Raster", "ANY"),
  function(
    formula, 
    data,  
    grid,
    covariates=NULL,
    buffer=0,
    ...) {

    dataCov = gm.dataRaster(
      formula, data,
      grid,
      covariates,
      buffer)

    callGeneric(
      formula = dataCov$formula, 
      data = dataCov$data, 
      grid = dataCov$grid, 
      covariates = dataCov$covariates, ...)
  }
  )


setMethod("glgm", 
  signature("formula", "Spatial", "Raster", "ANY"),
  function(formula, 
    data, grid, 
    covariates=NULL, 
    buffer=0,...) {

    dataCov = gm.dataSpatial(
      formula, data, 
      grid, covariates, buffer)
    callGeneric(formula, 
      data=dataCov$data@data, 
      grid=dataCov$grid, 
      covariates=dataCov$covariates, ...)
  }
  )

#################
#### the real work
##################

setMethod("glgm", 
  signature("formula", "data.frame", "Raster", "data.frame"), 
  function(formula, data,  grid, 
    covariates=NULL, 
    shape=1, priorCI=NULL, 
    mesh=FALSE,...) {

    if(!any(names(grid)=='space'))
      warning("grid must have a layer called space with inla cell ID's")

    allVars = allVarsP(formula)

    if(!all(allVars %in% names(data)) )
      warning("some covariates seem to be missing: formula ", 
        paste(allVars, collapse=" "), ", data: ", 
        paste(names(data), collapse=" "))

    cells = trim(grid[['space']])
    firstCell = values(cells)[1]
    cellDim = dim(cells)[1:2]
      # first cell = 2 * buffer^2 + ncolSmall * buffer + buffer
      # buffer = -(nrowSmall+1) + sqrt (  (nrowSmall+1)^2 + 8 firstCell / 4
    buffer = (-(cellDim[1]+1) + sqrt(  (cellDim[1]+1)^2 + 8* (firstCell-1) ))/4
      # data, cells, and covariates must have varilable called 'space'		
      # values of cells must be index numbers, and cells shouldnt include the buffer		
    forInla = thedots = list(...)

      # priors for spatial standard deviation and nugget std dev.
    sdNames = unique(c("sd",grep("^sd", names(priorCI), value=TRUE)))
      # if model is Gaussian, look for prior for sdNugget
    if(!any(names(thedots)=="family")) {
      forInla$family =  "gaussian"
    }
    if(thedots$family=="gaussian") {
      sdNames = unique(c(sdNames, "sdNugget"))
    }
    if(thedots$family=="gamma") {
      sdNames = unique(c(sdNames, "gammaShape"))
    }
    if(thedots$family %in% c("weibull", "weibullsurv") ) {
      sdNames = unique(c(sdNames, "weibullShape"))
    }

      # list of prior distributions
    if(any(names(priorCI)=='distributions')){
      priorDistributions = priorCI$distributions
    } else {
      priorDistributions = list()
    }

      # priors for sd's (and precisions) 

    precPrior=list()

    for(Dsd in sdNames) {
      Dprec = gsub("^sd","precision",Dsd)
      if(any(names(priorDistributions)==Dprec)) {
          # distribution supplied

        if('scale' %in% names(priorDistributions[[Dprec]])){
          priorDistributions[[Dprec]]['rate'] = 1/ priorDistributions[[Dprec]]['scale']
        }

        if(all(c('shape','rate') %in% names(priorDistributions[[Dprec]]))) {

          precPrior[[Dsd]] = list(
            params = c(
              shape=as.numeric(priorDistributions[[Dprec]]['shape']), 
              rate=as.numeric(priorDistributions[[Dprec]]['rate'])),
            prior = 'loggamma')

        } else {
          precPrior[[Dsd]] = list(params=c(
            shape=priorDistributions[[Dprec]][1],
            rate=priorDistributions[[Dprec]][2]),
          prior = 'loggamma')
        }
      } else if(any(names(priorCI)==Dsd)) {
          # if it's a matrix first column is sd, 2nd column is density
        if(is.matrix(priorCI[[Dsd]])) {
            # get rid of zero
          priorCI[[Dsd]] = priorCI[[Dsd]][
          priorCI[[Dsd]][,1]>0,]
          logPrecDens = cbind(
            -2 * log(priorCI[[Dsd]][,1]),
            priorCI[[Dsd]][,2] * priorCI[[Dsd]][,1]/2
            )
            # reorder smallest to largest precision		
          logPrecDens = logPrecDens[order(logPrecDens[,1]), ]

          precPrior[[Dsd]] = list(
            prior=paste("table:", 
              paste(as.vector(logPrecDens), collapse=' ')
              )
            )
          precPrior[[Dsd]]$string = paste(
            "prior='",
            precPrior[[Dsd]]$prior,
            "'", sep=''
            )
        } else { # not a table

        if(length(priorCI[[Dsd]])==1){
          priorCI[[Dsd]] = c(
            u = as.numeric(priorCI[[Dsd]]),
            alpha = 0.05
            )
        }
        if(priorCI[[Dsd]][2]<priorCI[[Dsd]][1]){
          names(priorCI[[Dsd]]) = c('u','alpha')
        }

            # find distribution from interval supplied
            # if of length 1, it's pc prior u with alpha = 0.05

        if(!length(names(priorCI[[Dsd]])))
          names(priorCI[[Dsd]]) = c('lower','upper')

        if(all(c('u','alpha') %in% names(priorCI[[Dsd]]))) {
              # pc priors
          precPrior[[Dsd]] = list(
            params=priorCI[[Dsd]],
            prior = 'pc.prec')
        } else {
              # gamma prior

          obj1 = sort(priorCI[[Dsd]]^-2)
          cifun = function(pars) {
            theci = 	pgamma(obj1, shape=pars[1], 
              rate=pars[2],log.p=T)

            (log(0.025) - theci[1])^2 +
            (2*(log(0.975) - theci[2]))^2		
          }
          precPrior2=optim(c(.5,.5/mean(obj1)), cifun, 
            lower=c(0.000001,0.0000001),method="L-BFGS-B")

          names(precPrior2$par) = c("shape","rate")

          precPrior[[Dsd]] = list(
            params = precPrior2$par,
            prior = 'loggamma')
        } # end gamma prior
      }
    } else { # no prior supplied
          # default prior
    precPrior[[Dsd]] = list(
      params = c(shape=0.01, rate=0.01),
      prior = 'loggamma')
  }
  if(!length(precPrior[[Dsd]]$string))
    precPrior[[Dsd]]$string = paste("param=c(",
      paste(precPrior[[Dsd]]$params, collapse=","),
      "), prior='", precPrior[[Dsd]]$prior, "'")
} # end loop Dsd



if(any(names(priorDistributions)=='range')) {
        # distribution supplied

  if('scale' %in% names(priorDistributions[['range']])){
    priorDistributions[['range']]['rate'] = 1/ priorDistributions[['range']]['scale']
  }


  if(all(c('shape','rate') %in% names(priorDistributions$range))) {

    ratePrior = c(
      shape=as.numeric(priorDistributions$range['shape']), 
      rate=as.numeric(priorDistributions$range['rate']*xres(cells))
      )

  } else {
    ratePrior = c(
      shape=priorDistributions$range[1],
      rate = priorDistributions$range[2]*xres(cells)
      )
  }
} else if("range" %in% names(priorCI)) {

  if(is.matrix(priorCI$range)) {

    priorCI$range = priorCI$range[
    priorCI$range[,1]>0,]

          # prior on log(1/(range/cellsize)) = log(cellsize) - log(range)
    logRangeDens = cbind(
      log(xres(cells))-log(priorCI$range[,1]),
      priorCI$range[,2] * priorCI$range[,1]
      )
          # reorder smallest to largest precision		
    logRangeDens = logRangeDens[nrow(logRangeDens):1, ]
    ratePrior = list(string=paste(
      "prior='table: ", 
      paste(as.vector(logRangeDens), collapse=' '),
      "'", sep=''
      ))
  } else if (priorCI$range[2] < priorCI$range[1]){ 
          # pc prior
    ratePrior = pcPriorRange(
      q = priorCI$range[1],
      p = priorCI$range[2],
      cellSize = cells
      )
  } else {
          # gamma prior
    if(priorCI$range[1] < xres(cells)/4) {
      priorCI$range[1] = xres(cells)/4
      warning("lower bound of range CI too small, setting it to 1/4 cell size")
    }

          # rang parameter, in terms of cells, not km.
    obj1=sort(priorCI$range/xres(cells))

    cifun = function(pars) {
      theci = 		pgamma(obj1, shape=pars[1], rate=pars[2], log.p=T)

      (theci[1] - log(0.025))^2 +
      (theci[2] - log(0.925))^2 
    }

    ratePrior2=optim(c(2,2/mean(obj1)), cifun, 
      lower=c(0.001,0.001),method="L-BFGS-B")
    ratePrior = ratePrior2$par
    names(ratePrior ) = c("shape","rate")
    ratePrior = list(
      params = ratePrior, prior='loggamma'
      )
  } 
}	else { # no range in priorCI
ratePrior = list(params=c(shape=0.01, rate=0.01), prior='loggamma')
}
if(!length(ratePrior$string))
  ratePrior$string = paste("param=c(",
    paste(ratePrior$params, collapse=","),
    "), prior='",ratePrior$prior,"'", sep='')

      # prior for gamma shape
      # log-normal, priorCI is 4 standard deviations
familyShapeName = grep("(familly|gamma|weibull)Shape", names(priorCI), value=TRUE)
if( length(familyShapeName) ) {
  familyShapeName = familyShapeName[1]
  familyShapePrior  = list(
    prior='gaussian',
    param=c(
      mean=as.numeric(mean(log(priorCI[[familyShapeName]]))),
      precision = as.numeric(abs(diff(log(priorCI[[familyShapeName]])))[1]/4)^(-2)
      )
    )
} else {
  familyShapePrior = NULL
}
      # done priors

      # formula for spatial random effect
spaceFormula = c(
  ".~.+ f(space, model='matern2d', ",
  "nrow=", nrow(cells)+2*buffer, 
  ", ncol=", ncol(cells)+2*buffer,
  ", nu=", shape, 
  ", hyper = list(",
  "range=list(",
  ratePrior$string,
  "),",
  "prec=list(",
  precPrior$sd$string,
  ")",
  " ) )", 
  sep=""
  )

formulaOrig = formula

if(any(names(forInla)=='spaceFormula')) {
  formula = update.formula(formula,	as.formula(forInla$spaceFormula))
  forInla = forInla[setdiff(names(forInla), 'spaceFormula')]
} else {
  formula = update.formula(formula,	as.formula(spaceFormula))
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
      covariates[,D] = factor(covariates[,D],
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
  theNA = apply(data[,allVars,drop=FALSE], 
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
if(!is.null(precPrior$sdNugget)) {
  forInla$control.family$hyper$prec =
  list(prior=precPrior$sdNugget$prior)
  if(length(precPrior$sdNugget$params))
    forInla$control.family$hyper$prec$params = 
  precPrior$sdNugget$params
}
if(length(familyShapeName)) {
  forInla$control.family$hyper$theta =
  familyShapePrior 
}


      # get rid of some elements of forInla that aren't required
forInla = forInla[grep("^buffer$", names(forInla), invert=TRUE)]
if(!length(forInla$lincomb)) forInla$lincomb = NULL


if(requireNamespace("INLA", quietly=TRUE)) {
  if(identical(forInla$verbose, TRUE)) {
    tFile = tempfile('glgm', tempdir(), '.rds')
    message(paste('saving INLA objects as', tFile))
    saveRDS(forInla, file=tFile)
  }
  inlaResult = try(do.call(INLA::inla, forInla))
} else {
  inlaResult = 
  list(logfile="INLA is not installed. \n install splines, numDeriv, Rgraphviz, graph,\n fields, rgl, mvtnorm, multicore, pixmap,\n splancs, orthopolynom \n then see www.r-inla.org")
}
if(identical(forInla$verbose, TRUE)) {
  message("inla done") 
}
if(all(names(inlaResult)=="logfile") | class(inlaResult) == 'try-error')
  return(c(forInla, inlares=inlaResult))


params = list(range=list())

      # posterior distributions forrange

params$range$posterior = 
inlaResult$marginals.hyperpar[["Range for space"]] %*% 
diag(c(xres(cells), 1/xres(cells)))

      # parameter priors for result
if(any(names(ratePrior) == 'priorRange')) {
        # pc prior
  params$range$params = ratePrior$lambda
  params$range$prior = ratePrior$priorRange
  params$range$priorScale = ratePrior$priorScale
  params$range$posteriorScale = cbind(
    x = 1/inlaResult$marginals.hyperpar[['Range for space']][,'x'],
    y = inlaResult$marginals.hyperpar[['Range for space']][,'y'] *
      inlaResult$marginals.hyperpar[['Range for space']][,'x']^(2)
    )
        # add third column to posterior matrix with prior
  params$range$posterior = cbind(
    x=params$range$posterior[,1],
    y=params$range$posterior[,2],
    prior = stats::dexp(
      xres(cells)/params$range$posterior[,1], 
      ratePrior$lambda) * (params$range$posterior[,1])^(-2) * 
    xres(cells)
    )

} else if(!is.matrix(priorCI$range)) {
  params$range$userPriorCI = priorCI$range
  params$range$priorCI = 
  xres(cells) *
  qgamma(c(0.025,0.975), 
    shape=ratePrior$params["shape"], 
    rate=ratePrior$params["rate"])
  params$range$priorCIcells = 
  qgamma(c(0.975,0.025), 
    shape=ratePrior$params["shape"], 
    rate=ratePrior$params["rate"])
  params$range$params.intern = ratePrior$params
  params$range$params = c(
    ratePrior$params['shape'],
    ratePrior$params['rate']/xres(cells))

  rangeLim = 	c(
    qgamma(c(0.001,0.999), 
      shape=params$range$params["shape"], 
      rate=params$range$params["rate"]),
    max(
      params$range$posterior[,1], na.rm=TRUE
      ))
  rangeSeq = c(0, seq(min(rangeLim), max(rangeLim), len=999))
  params$range$prior=cbind(
    x=rangeSeq,
    y=dgamma(rangeSeq,         
      shape=params$range$params["shape"], 
      rate=params$range$params["rate"])
    )
        # add third column to posterior matrix with prior
  params$range$posterior = cbind(
    x=params$range$posterior[,1],
    y=params$range$posterior[,2],
    prior = dgamma(
      params$range$posterior[,1],
      shape=params$range$params["shape"], 
      rate=params$range$params["rate"])
    )

  params$range$postDiv1000 = params$range$posterior %*%
  diag(c(1/1000, 1000, 1000))

} else {
  params$range$prior = priorCI$range
  params$range$postDiv1000 = params$range$posterior %*% 
  diag(c(1/1000, 1000))
}

for(Dsd in names(precPrior)) {
  if(is.matrix(priorCI[[Dsd]])) {
    params[[Dsd]] = list(
      prior = priorCI[[Dsd]]
      )
  } else if(precPrior[[Dsd]]$prior == 'loggamma'){
    params[[Dsd]] = list(userPriorCI=priorCI[[Dsd]], 
      priorCI = 1/sqrt(
        qgamma(c(0.975,0.025), 
          shape=precPrior[[Dsd]]$params["shape"], 
          rate=precPrior[[Dsd]]$params["rate"])),
      params.intern=precPrior[[Dsd]])

    precLim = 	qgamma(c(0.999,0.001), 
      shape=precPrior[[Dsd]]$params["shape"], 
      rate=precPrior[[Dsd]]$params["rate"])
    sdLim = 1/sqrt(precLim)
    sdSeq = seq(min(sdLim), max(sdLim), len=1000)
    precSeq = sdSeq^(-2)
    params[[Dsd]]$prior=cbind(
      x=sdSeq,
      y=dgamma(precSeq, shape=precPrior[[Dsd]]$params["shape"], 
        rate=precPrior[[Dsd]]$params["rate"]) *2* (precSeq)^(3/2) 
      )
  } else if( identical(precPrior[[Dsd]]$prior, 'pc.prec') ){
    params[[Dsd]] = list(userPriorCI=priorCI[[Dsd]], 
      priorCI = 1/sqrt(
        INLA::inla.pc.qprec(c(0.975,0.025),  
          u = precPrior[[Dsd]]$params['u'], 
          alpha = precPrior[[Dsd]]$params['alpha'])
        ),
      params.intern=precPrior[[Dsd]]$params)

    precLim = INLA::inla.pc.qprec(c(0.999,0.001),  
      u = precPrior[[Dsd]]$params['u'], 
      alpha = precPrior[[Dsd]]$params['alpha'])
    sdLim = 1/sqrt(precLim)
    sdSeq = seq(min(sdLim), max(sdLim), len=1000)
    precSeq = sdSeq^(-2)
    params[[Dsd]]$prior=cbind(
      x=sdSeq,
      y=INLA::inla.pc.dprec(precSeq, 
        u = precPrior[[Dsd]]$params['u'], 
        alpha = precPrior[[Dsd]]$params['alpha']
        ) * 2 * (precSeq)^(3/2) 
      )
  }
}

if(length(familyShapeName)) {

  paramsFamilyShape = 	c(
    familyShapePrior$param["mean"], 
    sd=as.numeric(1/sqrt(familyShapePrior$param["precision"]))
    )

  xLim = sort(exp(-stats::qnorm(
    c(0.999,0.001), 
    mean=paramsFamilyShape["mean"], 
    sd=paramsFamilyShape["sd"])
  ))

  xSeq  = c(0,seq(xLim[1], xLim[2], len=1000))

  params[[familyShapeName]] = list(
    userPriorCI = priorCI[[familyShapeName]],
    priorCI = sort(exp(-stats::qnorm(c(0.975,0.025), 
      mean=familyShapePrior$param["mean"], 
      sd=1/sqrt(familyShapePrior$param["precision"])))
    ),
    params.intern=familyShapePrior$param,
    params = paramsFamilyShape,
    distribution = 'lognormal',
    prior = cbind(
      x=xSeq, 
      y = stats::dlnorm(xSeq, meanlog = paramsFamilyShape['mean'],
        sdlog = paramsFamilyShape['sd'])
      )
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

  params[[Dsd]]$posterior=
  inlaResult$marginals.hyperpar[[thesd[Dsd]]]
  params[[Dsd]]$posterior[,"y"] = params[[Dsd]]$posterior[,"y"] * 2*  
  params[[Dsd]]$posterior[,"x"]^(3/2) 
  params[[Dsd]]$posterior[,"x"] = 1/sqrt(params[[Dsd]]$posterior[,"x"])  
  params[[Dsd]]$posterior = params[[Dsd]]$posterior[
  seq(dim(params[[Dsd]]$posterior)[1],1),]		

  if(identical(precPrior[[Dsd]]$prior, 'pc.prec') ) {

    params[[Dsd]]$posterior = cbind(
      x=params[[Dsd]]$posterior[,'x'],
      y=params[[Dsd]]$posterior[,'y'],
      prior = INLA::inla.pc.dprec(
        params[[Dsd]]$posterior[,'x']^(-2), 
        u = precPrior[[Dsd]]$params['u'], 
        alpha = precPrior[[Dsd]]$params['alpha']
        ) * 2 * (params[[Dsd]]$posterior[,'x']^(-2))^(3/2) 
      )
}


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


# put range in summary, in units of distance, not numbers of cells
thecolsFull =c("mean","sd",thecols,"mode") 
params$summary["range",thecolsFull]=				
xres(cells)*
inlaResult$summary.hyperpar[
"Range for space",
thecolsFull
]
if(params$summary["range", 'mean'] > 1000) {

  params$summary["range",thecolsFull]=				
  params$summary["range",thecolsFull] / 1000
  rownames(params$summary) = gsub(
    "^range$", "range/1000",
    rownames(params$summary)
    )
}
dimnames(params$summary) = lapply(dimnames(params$summary),
  function(qq) {
    qq=gsub("_", "\\\\textunderscore~", qq)
    qq=gsub("\\$", "\\\\textdollar~", qq)
    qq=gsub("<", "\\\\textless~", qq)
    qq=gsub(">", "\\\\textgreater~", qq)
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



resRaster=stack(resRasterRandom, resRasterFitted, cells)

result=list(inla=inlaResult,
  raster=resRaster,
  parameters=params
  )

result

}
)