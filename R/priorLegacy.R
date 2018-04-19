priorLegacy = function(priorCI, family, cellSize) {


 # priors for spatial standard deviation and nugget std dev.

	names(priorCI) = gsub("^sdNugget", "sdObs", names(priorCI))
	sdNames = unique(c("sd",grep("^sd", names(priorCI), value=TRUE)))
	if(family=="gaussian") {
		sdNames = unique(c(sdNames, "sdObs"))
	}


    # if model is Gaussian, look for prior for sdNugget
	if(family=="gamma") {
		sdNames = unique(c(sdNames, "gammaShape"))
	}
	if(family %in% c("weibull", "weibullsurv") ) {
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
					param = c(
						shape=as.numeric(priorDistributions[[Dprec]]['shape']), 
						rate=as.numeric(priorDistributions[[Dprec]]['rate'])),
					prior = 'loggamma')

			} else {
				precPrior[[Dsd]] = list(param=c(
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
				precPrior[[Dsd]]$string = paste0(
					"list(prior='",
					precPrior[[Dsd]]$prior,
					"')", sep=''
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

				expRate = -log(priorCI[[Dsd]]['alpha'])/priorCI[[Dsd]]['u']

				precPrior[[Dsd]] = list(
					param=priorCI[[Dsd]],
					prior = 'pc.prec',
					dprior = eval(parse(text=paste0('function(x) stats::dexp(x, rate=', expRate, ')'))),
					extra = list(mean = unname(1/expRate), rate=unname(expRate))
					)
				environment(precPrior[[Dsd]]$dprior) = baseenv()

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
					param = precPrior2$par,
					prior = 'loggamma')


				precPrior[[Dsd]]$extra = list(
					ciProb = c(0.025, 0.975),
					userPriorCI = priorCI[[Dsd]],
					priorCI = sort(1/sqrt(stats::qgamma(
						c(0.025, 0.975),
						shape=precPrior[[Dsd]]$param['shape'], 
						rate=precPrior[[Dsd]]$param['rate']))),
					optim = precPrior2
					)
#    densmat[, "y"] = densmat[, "y"] * 2 * densmat[, "x"]^(3/2)
#    densmat[, "x"] = 1/sqrt(densmat[, "x"])

				precPrior[[Dsd]]$dprior = eval(parse(text=paste0(
					'function(x) 2 * x^(-3) * stats::dgamma(x^(-2), shape=',
					precPrior[[Dsd]]$param['shape'],
					',rate=', precPrior[[Dsd]]$param['rate'], ')')))
				environment(precPrior[[Dsd]]$dprior) = baseenv()

			} # end gamma prior
		}
	} else { # no prior supplied
          # default prior
	precPrior[[Dsd]] = list(
		param = c(shape=0.01, rate=0.01),
		prior = 'loggamma')
	precPrior[[Dsd]]$dprior = eval(parse(text=paste0(
			'function(x) 2 * x^(3/2) * stats::dgamma(x^(-2), shape=',
			precPrior[[Dsd]]$param['shape'],
			',rate=', precPrior[[Dsd]]$param['rate'], ')')))
	environment(precPrior[[Dsd]]$dprior) = baseenv()


}
if(!length(precPrior[[Dsd]]$string))
	precPrior[[Dsd]]$string = paste0("list(param=c(",
		paste0(precPrior[[Dsd]]$param, collapse=","),
		"), prior='", precPrior[[Dsd]]$prior, "')")
} # end loop Dsd



if(any(names(priorDistributions)=='range')) {
        # distribution supplied

	if('scale' %in% names(priorDistributions[['range']])){
		priorDistributions[['range']]['rate'] = 1/ priorDistributions[['range']]['scale']
	}


	if(all(c('shape','rate') %in% names(priorDistributions$range))) {

		ratePrior = c(
			shape=as.numeric(priorDistributions$range['shape']), 
			rate=as.numeric(priorDistributions$range['rate']*cellSize)
			)

	} else {
		ratePrior = c(
			shape=priorDistributions$range[1],
			rate = priorDistributions$range[2]*cellSize
			)
	}
} else if("range" %in% names(priorCI)) {

	if(is.matrix(priorCI$range)) {

		priorCI$range = priorCI$range[
		priorCI$range[,1]>0,]

          # prior on log(1/(range/cellsize)) = log(cellsize) - log(range)
		logRangeDens = cbind(
			log(cellSize)-log(priorCI$range[,1]),
			priorCI$range[,2] * priorCI$range[,1]
			)
          # reorder smallest to largest precision		
		logRangeDens = logRangeDens[nrow(logRangeDens):1, ]
		ratePrior = list(string=paste0(
			"list(prior='table: ", 
			paste0(as.vector(logRangeDens), collapse=' '),
			"')", sep=''
			))

		ratePrior$dprior = list(
			range = stats::approxfun(priorCI$range[,1], priorCI$range[,2]),
			scale = stats::approxfun(1/priorCI$range[,1], priorCI$range[,1]^(-2)*priorCI$range[,2]) 
			)

	} else if (priorCI$range[2] < priorCI$range[1]){ 
          # pc prior
		ratePrior = pcPriorRange(
			q = priorCI$range[1],
			p = priorCI$range[2],
			cellSize = cellSize
			)
	} else {
          # gamma prior
		if(priorCI$range[1] < cellSize/4) {
			priorCI$range[1] = cellSize/4
			warning("lower bound of range CI too small, setting it to 1/4 cell size")
		}

          # rang parameter, in terms of cells, not km.
		obj1=sort(priorCI$range/cellSize)

		cifun = function(pars) {
			theci = stats::pgamma(obj1, shape=pars[1], rate=pars[2], log.p=T)

			(theci[1] - log(0.025))^2 +
			(theci[2] - log(0.925))^2 
		}

		ratePrior2=optim(c(2,2/mean(obj1)), cifun, 
			lower=c(0.001,0.001),method="L-BFGS-B")
		ratePrior = ratePrior2$par
		names(ratePrior ) = c("shape","rate")
		ratePrior = list(
			param = ratePrior, prior='loggamma'
			)

		ratePrior$extra = list(
			userParam = ratePrior$param[c('shape','rate')] / c(1,cellSize),
			ciProb = c(0.025, 0.975),
			userPriorCI = priorCI$range,
			priorCI = cellSize * stats::qgamma(
				c(0.025, 0.975),
				shape=ratePrior$param['shape'], 
				rate=ratePrior$param['rate']),
			optim = ratePrior2
			)

		ratePrior$dprior = list(
			scale = eval(parse(text=paste0(
				'function(x) x^(-2)*stats::dgamma(1/x, shape=',
				ratePrior$extra$userParam['shape'],',rate=', 
				ratePrior$extra$userParam['rate'], ')'))),
			range = eval(parse(text=paste0(
				'function(x) stats::dgamma(x, shape=',
				ratePrior$extra$userParam['shape'],',rate=', 
				ratePrior$extra$userParam['rate'], ')')))
			)
		environment(ratePrior$dprior$range) = baseenv()
		environment(ratePrior$dprior$scale) = baseenv()
	} # end types of rate prior
} # end distribution for range not supplied

if(!length(ratePrior$string))
	ratePrior$string = paste0("list(param=c(",
		paste0(ratePrior$param, collapse=","),
		"), prior='",ratePrior$prior,"')", sep='')

result = c(
	precPrior,
	list(range = ratePrior)
	)


      # prior for gamma shape
      # log-normal, priorCI is 4 standard deviations

familyShapeName = grep("(family|gamma|weibull)Shape", names(priorCI), value=TRUE)
if( length(familyShapeName) ) {
	familyShapeName = familyShapeName[1]
	familyShapePrior  = list(
		prior='gaussian',
		param=c(
			mean=as.numeric(mean(log(priorCI[[familyShapeName]]))),
			precision = as.numeric(abs(diff(log(priorCI[[familyShapeName]])))[1]/4)^(-2)
			)
		)
	result$family =familyShapePrior
} 


result
}