spatialRocPolyTemplate = function(
		truth, fit
) {
	

	toRasterize = fit[[1]]$data
	toRasterize$fitID = 1:length(toRasterize)
	template = rasterize(
			toRasterize,
  		truth,
			field='fitID'
			)

	names(template) = 'fitID'
	
	template
	
}

spatialRocRasterTemplate = function(
		truth, fit
) {
	if(length(grep("^Raster", class(fit)))) {
		template = raster(fit)
	} else {
		template = raster(fit[[1]]$raster)
	}
	values(template) = seq(1, ncell(template))
	names(template) = 'fitID'
	
	# remove cells with no predictions

	if(length(grep("^Raster", class(fit)))) {
		template = mask(template, fit[[1]])
	} else {
		template = mask(template, fit[[1]]$raster$predict.mean)
	}
	
	templateID = stackRasterList(
			list(fitID = template), truth, method='ngb'
	)

	# and cells with no truth
	templateID = mask(templateID, truth[[1]])
	
	templateID
}

spatialRocSims = function(
		truthCut, marginals, templateID, 
		breaks, prob
) {
	
	
	breaksInterior = breaks[seq(2, length(breaks)-1)]
	
	
	Slevels = seq(from=1.5,by=1,len=length(breaks)-1)
	
	SlevelsRound = 1+1:length(Slevels)
	names(Slevels) = SlevelsRound
	
	# breaks = -Inf, 0.4, 0.7, 0.8, Inf
	# breaksInterior = 0.4, 0.7, 0.8
	# SlevelsRound = 1 2 3 4
  # Slevels = 1.5, 2.5, 3.5, 4.5
  # There are 4 bins, 5 breaks, 3 interior breaks
	# Slevels are bin centres
	# SlevelsRound are bin numbers
	
	
	if(is.null(names(marginals)))
		names(marginals) = 1:length(marginals)
		
	templateID = mask(templateID, truthCut[[1]])
	
	truthCdf = zonal(
			truthCut,
			templateID,
			function(x,...) {
				stats::ecdf(x)(Slevels)
			}
	)
	
	
	rownames(truthCdf) = truthCdf[,'zone']
	truthCdf = truthCdf[,grep('zone', colnames(truthCdf), invert=TRUE)]
	colnames(truthCdf) = paste(
			rep(names(marginals), 
					rep(length(Slevels),nlayers(truthCut))
			), 
			rep(Slevels, nlayers(truthCut)),
			sep="_"
	)
	truthCdf = array(
			truthCdf,
			c(nrow(truthCdf), length(Slevels), nlayers(truthCut)),
			dimnames = list(
					rownames(truthCdf),
					Slevels,
					names(marginals)
			)
	)
	truthCdfUpper = 1-truthCdf
	
	idTable = table(values(templateID))		
	idMatrix = matrix(
			idTable[dimnames(truthCdf)[[1]]],
			length(idTable),
			length(prob)
	)
	
	truthCells = array(
			idTable[dimnames(truthCdf)[[1]]],
			dim=dim(truthCdf)
	)
	
	truthOver = truthCdfUpper * truthCells 
	truthUnder = truthCdf * truthCells 
	
	
	
	allP = allN = tP = tN = fP = fN = array(NA,
			c(length(prob), length(SlevelsRound), length(marginals)),
			dimnames= list(
					prob,SlevelsRound, names(marginals)
			)
	)
	
		for(Dsim in names(marginals)) {
		
		

		if(length(grep("^Raster", class(marginals[[Dsim]])))) {		
			# local-em
			
			pMat = cbind('1'=1,as.data.frame(
							marginals[[Dsim]][[grep("threshold.", names(marginals[[Dsim]]))]]
					), '0'=0)
			colnames(pMat) = c('1',names(Slevels))
			pMat = pMat[sort(na.omit(unique(values(templateID)))),]
		} else { # not local-em

		notNull = dimnames(truthCdf)[[1]]
		if(! table(notNull %in% names(marginals[[Dsim]])) )
			warning("can't map prediction ID's to marginal distributions")

			# pMat is prob below break

	  	# don't include last break (infinity) or first break
			pMat = simplify2array(
					lapply(
							marginals[[Dsim]][notNull], 
							INLA::inla.pmarginal, 
							q=breaksInterior)
			)
			if(is.vector(pMat))
				pMat = t(as.matrix(pMat))
			rownames(pMat) = seq(2, by=1, len=nrow(pMat))
			# upper tail probabilities
			pMat = pMat[,match(names(marginals[[Dsim]]), colnames(pMat)), drop=FALSE]
			pMat = t(pMat)
			
			pMat = pMat[as.numeric(dimnames(truthCdf)[[1]]),,drop=FALSE]
			# pMat is prob above break
			pMat = 1-pMat
			pMat = cbind('1'=1, pMat,0)
			colnames(pMat)[ncol(pMat)] = as.character(ncol(pMat))
		}
	
		for(Dbin in names(Slevels)) {
			
			Dmidpoint = as.character(Slevels[Dbin])
			
			predOver = outer(pMat[,Dbin], prob, '>')
#					colnames(predOver) = prob
			
			tpMat = predOver*truthOver[,Dmidpoint,Dsim]
			fpMat = predOver*truthUnder[,Dmidpoint,Dsim]
			tnMat = (1-predOver)*truthUnder[,Dmidpoint,Dsim]
			fnMat = (1-predOver)*truthOver[,Dmidpoint,Dsim]
			allnMat = tnMat + fpMat
			allpMat = tpMat + fnMat
			
			tP[,Dbin,Dsim] = apply(tpMat,2,sum, na.rm=TRUE)
			tN[,Dbin,Dsim] = apply(tnMat,2,sum, na.rm=TRUE)
			fP[,Dbin,Dsim] = apply(fpMat,2,sum, na.rm=TRUE)
			fN[,Dbin,Dsim] = apply(fnMat,2,sum, na.rm=TRUE)
			allP[,Dbin,Dsim] = apply(allpMat,2,sum, na.rm=TRUE)
			allN[,Dbin,Dsim] = apply(allnMat,2,sum, na.rm=TRUE)
			
		}
		
		
	}
	
	
	result = list(
			tP = tP, 
			allP = allP, 
			tN = tN, 
			allN = allN
	)
	result$allP[result$allP==0] = NA
	result$allN[result$allN==0] = NA
	
	return(result)
	
}


spatialRoc = function(fit, 
		rr=c(1,1.2, 1.5,2), truth, 
		border=NULL, random=FALSE,
		prob = NULL, spec = seq(0,1,by=0.01)){
	
	if(is.null(prob)){
		prob = 2^seq(-1, -10, len=20)
	}
	prob = sort(unique(c(prob, 1-prob, 0.5)))
	
	if(any(names(fit)=='inla') | length(grep("^Raster", class(fit)))){
		fit = list(fit)
	}

	thresholdNames = grep("^threshold\\.[[:digit:]]", names(fit[[1]]), value=TRUE)
	if(length(thresholdNames)) {
		# this is a local em result
		# rr must be that used for the bootstrap simulations
		rr = as.numeric(gsub("^threshold\\.", "", thresholdNames))
		random = FALSE
		isLocalEm = TRUE
	} else {
		isLocalEm = FALSE
	}
	breaks = c(-Inf, log(rr), Inf)
	
	
	if('raster'%in% names(truth))
		truth = truth$raster
	
	if(random) {
		truthVariable = 'random'
		truth = truth[[
				intersect(
						paste(truthVariable, c('',1:length(fit)), sep=''),
						names(truth)
				)
		]]
	} else { 
		truthVariable = 'relativeIntensity'
		
		truth = truth[[
				intersect(
						paste(truthVariable, c('',1:length(fit)), sep=''),
						names(truth)
				)
		]]
		tname = names(truth)
		truth = log(truth)
		names(truth) = tname
	} 

	if(is.null(border))
		if(any(names(fit[[1]])=='data'))
			border = fit[[1]]$data
		
	if(!is.null(border))
		truth = mask(truth, border)

	truthCut = cut(truth, breaks=breaks)
	names(truthCut) = names(truth)
	
	
	if(isLocalEm) {
		templateID = spatialRocRasterTemplate(
				truthCut, fit[[1]]
		)
		marginals = fit
	 		
	} else if(any(names(fit[[1]]) =='raster')){
		# lgcp or glgm
		
		templateID = spatialRocRasterTemplate(
			truthCut, fit
		) 
		
		
		if(random) {
			marginals = lapply(
				fit, function(x){
					 x=x$inla$marginals.random$space
					 names(x) = 1:length(x)
					 x
				}
				)
		} else {
			marginals = lapply(
				fit, function(x){
					x=x$inla$marginals.predict
					names(x) = 1:length(x)
					x
				}
			)
		}
		
	} else { # bym model

		templateID = spatialRocPolyTemplate(
				truthCut, fit
		) 
		
		if(random) {
			marginals = lapply(
					fit, function(x){
					  x$inla$marginals.bym
					}
			)
		} else {
			marginals = lapply(
					fit, function(x){
					  x$inla$marginals.fitted.bym
					}
			)
		}
		
	}

	res = spatialRocSims(
			truthCut, marginals, templateID, 
			breaks, prob
		)	
	
	

			bySim = abind::abind(
			onemspec = 1-res$tN / res$allN,
			sens = res$tP / res$allP,
			along=4
			)
			bySim = bySim[,-dim(bySim)[2],,,drop=FALSE]

			result = apply(bySim, c(1,2,4), mean)

			if(length(rr)>1) {
				dimnames(result)[[2]] = rr
			} 
			
			if(!is.null(spec) & length(rr)>1) {
				
				resultOut = matrix(
						NA,
						length(spec), 
						dim(result)[2]+1,
						dimnames = list(
								1:length(spec),
								c('onemspec', dimnames(result)[[2]])
								)
						)
				resultOut[,'onemspec'] = 1-spec
				for(D in dimnames(result)[[2]]) {
				  
				  resultOneMSpec = result[,D,'onemspec']
				  resultSens = result[,D,'sens']
				  
				  if(all(is.na(resultOneMSpec)) | 
				      all(is.na(resultSens)) | 
				      length(unique(resultOneMSpec[!is.na(resultOneMSpec)])) == 1) {
				    resultOut[,D] = NA
				  } else {
				    resultOut[,D] = approx(
				      resultOneMSpec, 
				      resultSens, 
				      xout = resultOut[,'onemspec'])$y
				  }
				}
				resultOrig = result
				result = resultOut
				attributes(result)$orig = resultOrig
			} else {
				result = drop(result)
			}
			
			attributes(result)$sim = bySim
	
  result
}
