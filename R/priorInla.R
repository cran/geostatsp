precToSd = function (densmat) 
{
	densmat[, "y"] = densmat[, "y"] * 2 * densmat[, "x"]^(3/2)
	densmat[, "x"] = 1/sqrt(densmat[, "x"])

	theMaxX = max(densmat[densmat[,'y'] > 10^(-5) * max(densmat[,'y']),'x'])
	do.call(cbind, approx(
		densmat[,'x'],   
		densmat[,'y'], 
		seq(0, theMaxX, len=nrow(densmat)),
		rule = 2
		))
}

priorInla = function(x, family='gaussian', cellSize=1) {

	names(x) = gsub("^sdNugget$", "sdObs", names(x))

	if("range" %in% names(x)) {
			# default prior if none specified
			# pc prior with median 10 cells
		if(!length(x$range)) 
			x$range = c(u=10*cellSize, alpha=0.5)
	}

	for(D in names(x)) {
		if(!is.list(x[[D]])) {
			if(length(x[[D]])==1) {
				# pc prior, median is specified
				x[[D]] = c(u=unname(x[[D]]), alpha = 0.5)
			}
			if(is.null(names(x[[D]])) & !is.matrix(x[[D]])) {
				if(x[[D]][2] > 1) {
					names(x[[D]]) = c('lower', 'upper')
				} else {
					names(x[[D]]) = c('u', 'alpha')
				}
			}
		}
	}

	sdNames = unique(c("sd",grep("^sd", names(x), value=TRUE)))
	if(family=="gaussian") {
		sdNames = unique(c(sdNames, "sdObs"))
	}

    # if model is Gaussian, look for prior for sdNugget
	if(family=="gamma") {
		sdNames = unique(c(sdNames, "gammaShape"))
	}

	ciTarget = log(c(0.025, 0.975))
	parGetRid = c('cellSize','info','dprior', 'extra')


	precPrior=list()

	for(Dsd in sdNames) {

		if(is.null(x[[Dsd]])) { 
		# no prior supplied
          # default prior
			x[[Dsd]] = c(u=1, alpha = 0.5)
		} 

		if(is.list(x[[Dsd]])) { 
	# list provided, use it unchanged
			precPrior[[Dsd]] = x[[Dsd]]
			precPrior[[Dsd]]$info = 'user supplied prior'
			if(!is.null(precPrior[[Dsd]]$initial))
				precPrior[[Dsd]]$initial = precPrior[[Dsd]]$initial^(-2)
		} else if(is.matrix(x[[Dsd]])) {
			# matrix with numerical values of density provided
  #   first column is sd, 2nd column is density
			precPrior[[Dsd]] = list(
				dprior = stats::approxfun(x[[Dsd]][,1], x[[Dsd]][,2]))

  # get rid of zero
			x[[Dsd]] = x[[Dsd]][x[[Dsd]][,1]>0,]
			logPrecDens = cbind(
				-2 * log(x[[Dsd]][,1]),
				x[[Dsd]][,2] * x[[Dsd]][,1]/2
				)
  # reorder smallest to largest precision		
			logPrecDens = logPrecDens[order(logPrecDens[,1]), ]

#			precPrior[[Dsd]] = x[[Dsd]]
			precPrior[[Dsd]]$string = paste0(
				"list(prior='table:", 
				paste(as.vector(logPrecDens), collapse=' '), 
				"')", sep='')
		} else if(all(c('u','alpha') %in% names(x[[Dsd]]))) {

              # pc prior specified

			expRate = -log(x[[Dsd]]['alpha'])/x[[Dsd]]['u']
			precPrior[[Dsd]] = list(
				param = x[[Dsd]],
				info = 'exponential prior on standard deviation',
				prior = 'pc.prec',
				dprior = eval(parse(text=paste0('function(x) stats::dexp(x, rate=', expRate, ')'))),
				extra = list(mean = unname(1/expRate), rate=unname(expRate))
				)
			environment(precPrior[[Dsd]]$dprior) = baseenv()
		} else if(all(c('lower','upper') %in% names(x[[Dsd]]))) {
# quantiles provided
              # gamma prior for sd
			cifun = function(pars) {
				theci = pgamma(
					x[[Dsd]],
					shape=pars[1], 
					rate=pars[2],log.p=TRUE)

				sum(c(1,4)*(ciTarget - theci)^2)
			}

			precPrior2=stats::optim(
				c(.5,.5/sqrt(prod(x[[Dsd]]))), 
				cifun, 
				lower=c(0.000001,0.0000001),method="L-BFGS-B")
			names(precPrior2$par) = c("shape","rate")


			precPrior[[Dsd]] = gammaSd(precPrior2$par)

			precPrior[[Dsd]]$extra = list(
				ciProb = exp(ciTarget),
				userPriorCI = x[[Dsd]],
				priorCI = stats::qgamma(
					exp(ciTarget),
					shape=precPrior[[Dsd]]$param['shape'], 
					rate=precPrior[[Dsd]]$param['rate']),
				optim = precPrior2
				)

		} # end interval supplied

		# create inla string if it doesn't yet exist
		if(!length(precPrior[[Dsd]]$string)) {
			precPrior[[Dsd]]$string = deparse(
				precPrior[[Dsd]][setdiff(names(precPrior[[Dsd]]), parGetRid)], 
				control=NULL)
		}

	} # end loop Dsd

# range
	rangePrior = NULL
	if("range" %in% names(x)) {
			# default prior if none specified
			# pc prior with median 10 cells

		if(is.list(x$range)) {
			if(is.null(x$range$cellSize))
				x$range$cellSize = cellSize
			if(!is.null(x$range$initial))
				x$range$initial = x$range$cellSize/x$range$initial
			if(identical(x$range$prior, 'invgamma')){
				rangePrior = gammaScale(x$range$param, cellSize)
			} else {
			rangePrior = x$range
			rangePrior$dprior = list(
				range = function(x) {rep(NA, length(x))}
				)
			environment(rangePrior$dprior$range) = baseenv()
			rangePrior$dprior$scale=rangePrior$dprior$range
			}
		}

		if(is.matrix(x$range)) {
				# table supplied

			rangePrior = list(
				dprior = list(
					range = stats::approxfun(x$range[,1], x$range[,2]),
					scale = stats::approxfun(1/x$range[,1], x$range[,1]^(-2)*x$range[,2]) ),
				cellSize = cellSize,
				info = 'table provided for range')

			x$range = x$range[x$range[,1]>0,]

          # prior on log(1/(range/cellsize)) = log(cellsize) - log(range)
			logRangeDens = cbind(
				log(cellSize)-log(x$range[,1]),
				x$range[,2] * x$range[,1])

          # reorder smallest to largest
			logRangeDens = logRangeDens[nrow(logRangeDens):1, ]
			rangePrior$string = paste0(
				"list(prior='table: ", 
				paste0(as.vector(logRangeDens), collapse=' '),
				"')", sep='')

		} else if(all(c('u', 'alpha') %in% names(x$range) ) ) { 

    # pc prior specified
			rangePrior = pcPriorRange(
				q = x$range['u'],
				p = x$range['alpha'],
				cellSize = cellSize)

		} else if(!is.list(x$range)){
			rangePrior = list()
		        # gamma prior for scale
			scaleTarget = unname(sort(cellSize/x$range[c('lower','upper')]))
			cifun = function(pars) {
				theci = stats::pgamma(
					scaleTarget,
					shape=pars[1], 
					rate=pars[2], log.p=TRUE)

				sum(c(1,4)*(ciTarget - theci)^2)
			}

			scalePrior2=stats::optim(
				c(.5,.5/sqrt(prod(scaleTarget))), 
				cifun, 
				lower=c(0.000001,0.0000001),method="L-BFGS-B")

			names(scalePrior2$par) = c("shape","rate")

			rangePrior = gammaScale(
				scalePrior2$par, cellSize
				)


			rangePrior$extra = c(
				rangePrior$extra,
				list(
				ciProb = exp(ciTarget),
				userPriorCI = x$range[c('lower','upper')],
				priorCI = sort(cellSize/stats::qgamma(
					exp(ciTarget),
					shape=rangePrior$param['shape'], 
					rate=rangePrior$param['rate'])),
				optim = scalePrior2
				))

		} # end types of prior

		if(!length(rangePrior$string))
			rangePrior$string = deparse(rangePrior[setdiff(names(rangePrior), parGetRid)], control=NULL)

	} # end have range

      # prior for gamma shape
      # log-normal, priorCI is 4 standard deviations

	result = c(
		precPrior,
		list(range = rangePrior)
		)


	familyShapeName = grep(
		"(family|gamma|weibull)Shape", names(x), value=TRUE)


	if( length(familyShapeName) ) {
		familyShapeName = familyShapeName[1]
		if(is.list(x[[familyShapeName]])) {
			familyShapePrior = list(
				string = deparse(x[[familyShapeName]][setdiff(names(x[[familyShapeName]]), parGetRid)]),
				name = familyShapeName
				)

			if(identical(x[[familyShapeName]]$prior, 'loggamma')) {
				familyShapePrior$dprior = eval(parse(text=paste(
					'function(x) dgamma(x, shape = ',
					x[[familyShapeName]]$param[1],
					', rate = ',
					x[[familyShapeName]]$param[2],					
					')') ))
			} 
		} else {
			# default prior
			familyShapePrior  = list(
				prior='gaussian',
				param=c(
					mean=as.numeric(mean(log(x[[familyShapeName]]))),
					precision = as.numeric(abs(diff(log(x[[familyShapeName]])))[1]/4)^(-2)
					),
				name = familyShapeName
				)
	}
		result$family = familyShapePrior
	}

	result
}
