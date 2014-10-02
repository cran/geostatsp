# the function to get gradients
grFun = function(parname, right=TRUE, 
		allParam, paramsRight, paramsLeft, observations,
		covMatOrig, covariates, obsBC,twoLogJacobian,
		coordinates, fixedVariance,reml,
		minustwotimes) {
	
	if(right) {
		allParam[parname] = paramsRight[parname]		
	} else {
		allParam[parname] = paramsLeft[parname]
	}
	
	if(parname=="boxcox") {
		obsBC =	((observations^allParam["boxcox"]) - 1)/
				allParam["boxcox"]
		twoLogJacobian = 2*(allParam["boxcox"]-1)* 
				sum(log(observations))					
	} else if(parname=="nugget") {
		covMat =covMatOrig
		Matrix::diag(covMat) = Matrix::diag(covMat) + allParam["nugget"]
	} else {
		covMat = matern(x=coordinates, param=allParam)	
		Matrix::diag(covMat) = Matrix::diag(covMat) + allParam["nugget"]
	}
	
	cholCovMat = Matrix::chol(covMat)
	cholCovInvX = Matrix::solve(cholCovMat, covariates)
	cholCovInvXcross = Matrix::crossprod(cholCovInvX)
	cholCovInvXcrossInv = Matrix::solve(cholCovInvXcross)
	cholCovInvY = Matrix::solve(cholCovMat, obsBC)
	betaHat = cholCovInvXcrossInv %*% Matrix::crossprod(cholCovInvX, cholCovInvY) 
	resids = obsBC - covariates %*% betaHat
	cholCovInvResid = Matrix::solve(cholCovMat, resids)
	totalSsq = as.vector(Matrix::crossprod(cholCovInvResid))
	totalVarHat = totalSsq/length(observations)
	
	if(!fixedVariance) { # profile likelihood with optimal sigma
		
		minusTwoLogLik = length(observations) * log(2*pi) + 
				length(observations) * log(totalVarHat) +
				2*sum(log(Matrix::diag(cholCovMat))) +
				length(observations) - twoLogJacobian
		
		if( reml ) {
			# REML
			choldet =  2*sum(log(Matrix::diag(chol(cholCovInvXcross)))) 
			
			minusTwoLogLik =  minusTwoLogLik + choldet
			
		}  # end reml
		
	} else { # a variance was supplied
		
		minusTwoLogLik = length(observations) * log(2*pi) +
				2*sum(log(Matrix::diag(cholCovMat))) +
				totalSsq - twoLogJacobian
		totalVarHat = 1
		if( reml ) {
			reml=FALSE
			warning("a variance was supplied by REML was requested.  doing ML instead.")		
		}
	}
	
	result = minusTwoLogLik
	if(minustwotimes) {
		names(result) = "minusTwoLogLik"
	} else {
		result = -0.5*result
		
		names(result)="logLik"
	} 
	if(reml)
		names(result) = gsub("Lik$", "RestrictedLik", names(result))
	result
}



logLikGradient = function(param, 
		observations, covariates, coordinates,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=NULL,
		lowerG = rep(-Inf, length(param)), 
		upperG=rep(Inf, length(param)), 
		eps=rep(0.001, length(param)),
		mc.cores=getOption("mc.cores", 1L)
	){
		upper = upperG
		lower = lowerG
		
		param=param[!is.na(param)]
		
		if(any(names(param)=="variance")) {
			warning("ignoring variance component of param, it will be estimated analytically")
			param = param[!names(param) %in% "variance"]
		}

		estNugget = any(names(param)=="nugget")
		if(any(names(moreParams)=="nugget")) {
			if(moreParams["nugget"]==0) {
				moreParams = moreParams[!names(moreParams)=="nugget"]		
			}
		}
		haveNugget = estNugget | any(names(moreParams)=="nugget")

		# is the variance parameter fixed?
		fixedVariance = any(names(moreParams)=="variance")	

# store box-cox if necessary

	estBoxCox = any(names(param)=="boxcox")
	if(any(names(moreParams)=="boxcox")) { 
		if(abs(moreParams["boxcox"]-1) < 0.001) {
			moreParams = moreParams[!names(moreParams)=="boxcox"]		
		}
	}
			
	useBoxCox = estBoxCox | any(names(moreParams)=="boxcox")
	if(useBoxCox) {
		if(!estBoxCox) { # box cox is fixed, transformed data must be passed in.
			if(any(names(stored)=="boxcox") ) {
				if(abs(stored$boxcox-moreParams["boxcox"])>0.0001)
					warning("boxcox param and stored are ", moreParams["boxcox"],
						" and ", stored$boxcox)

				twoLogJacobian = stored$twoLogJacobian
				obsBC = stored$observations
			} else { # compute stuff
				if(abs(moreParams["boxcox"]> 0.001)) {
					obsBC =	((observations^moreParams["boxcox"]) - 1)/
							moreParams["boxcox"]
				} else {
					moreParams["boxcox"] = 0
					obsBC = log(observations)					
				}
				
				twoLogJacobian = 2*(moreParams["boxcox"]-1)* 
						sum(log(observations))					
			}
		} else { # box cox is estimated, must compute
			if(abs(param["boxcox"]) > 0.00001) { 
			obsBC =	((observations^param["boxcox"]) - 1)/
							param["boxcox"]
				} else {
					obsBC = log(observations)
				}
			twoLogJacobian = 2*(param["boxcox"]-1)* 
					sum(log(observations))					
		}
	} else { # not useing box cox
		obsBC = observations
		twoLogJacobian = 0
	}		
		
	allParam = c(param, moreParams[!names(moreParams) %in% names(param)])
	
	# variance matrix
	covMatOrig = covMat = matern(x=coordinates, param=allParam)	
	if(haveNugget)
		Matrix::diag(covMat) = Matrix::diag(covMat) + allParam["nugget"]

	paramsRight = pmin(param + eps, upper[names(param)])
	paramsLeft = pmax(param - eps, lower[names(param)])
	



	# calcluate the gradient
	parForMapply = rep(names(param), 2)
	rightForMapply = rep(c(TRUE, FALSE), rep(length(param),2))

	if(length(observations)>200 ) {
		ncores = mc.cores
	} else {
		ncores = 1		
	}
	cat(ncores)

	moreArgs = list(obsBC=obsBC,twoLogJacobian=twoLogJacobian,
	allParam=allParam, paramsRight=paramsRight, 
	paramsLeft=paramsLeft, observations=observations,
	covMatOrig=covMatOrig,
	coordinates=coordinates, fixedVariance=fixedVariance,
	reml=reml, covariates=covariates,
	minustwotimes=minustwotimes)
	
if(ncores>1 & requireNamespace("parallel", quietly=TRUE)) {
	thegrad = parallel::mcmapply(grFun, parname=parForMapply, right=rightForMapply,
			MoreArgs=moreArgs,
			mc.cores=ncores)
} else {
	thegrad = mapply(grFun, parname=parForMapply, right=rightForMapply,
			MoreArgs=moreArgs)
}
	
	thegrad = matrix(thegrad, nrow=length(names(param)), ncol=2, 
					dimnames = list(names(param), 
							c("left","right")[1+rightForMapply[c(1,length(param)+1)]]))
			
	thegrad = thegrad[,"right"] - thegrad[,"left"]
	thegrad = thegrad * (paramsRight[names(thegrad)] - paramsLeft[names(thegrad)])
	thegrad
	
}

loglikLgmG = function(param,
		observations, covariates, coordinates,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=NULL,
		lowerG=NULL, upperG=NULL, 
		eps=NULL) {

	loglikLgm(param,
			data=observations, 
			formula=covariates, coordinates=coordinates,
			reml=reml, 
			minustwotimes= minustwotimes,
			stored=stored, moreParams=moreParams) 
	
}	
	

likfitLgmG = function(
		data, trend, 
		coordinates=data,
		param=c(range=1,nugget=0,shape=1),
		upper=NULL,lower=NULL, parscale=NULL,
		paramToEstimate = c("range","nugget"),
		reml=TRUE) {
	
	# for some reason thing break if I remove this next line...
	stuff = (class(coordinates))
	
	# check for the variance parameter
	estimateVariance = TRUE
	if(any(paramToEstimate=="variance")) {
		# remove varinace, it's estimated by profile likelihood
		paramToEstimate = paramToEstimate[paramToEstimate != "variance"]
		param = param[names(param)!="variance"]
	} else {
		if(any(names(param)=="variance")){
			estimateVariance = FALSE
			warning("variance is fixed and not estimated. If this isn't what you wanted remove variance from param")
		}		
	}
	
	# limits
	lowerDefaults = c(nugget=0,range=0,anisoRatio=0.0001,
			anisoAngleRadians=-pi/2,anisoAngleDegrees=-90,
			shape=0.01,boxcox=-3,variance=0)
	
	upperDefaults= c(nugget=Inf,range=Inf,anisoRatio=Inf,
			anisoAngleRadians=pi/2,anisoAngleDegrees=90,
			shape=0,boxcox=3,variance=Inf)
	
	
	lowerDefaults[names(lower)]=lower
	
	upperDefaults[names(upper)] = upper
	
	# par scale
	parscaleDefaults = c(range=NA, #range is delt with below
			nugget=1,
			boxcox=1,
			anisoAngleDegrees=10,
			anisoAngleRadians=2,
			anisoRatio=1,
			variance=1)
	
#	if(length(grep("^SpatialPoints", class(coordinates))))
#		print("wah")
	
	# convert input data to a model matrix
	if(class(trend)=="formula") {
		
		data = as.data.frame(data)
		theNA = apply(
				data[,all.vars(formulaRhs(trend)),drop=FALSE],
				1, function(qq) any(is.na(qq)))
		noNA = !theNA
		
		covariates = model.matrix(trend, data[noNA,])
		observations = formulaLhs(trend)
		
		if(!any(names(data)==observations))
			warning("can't find observations ", observations, "in data")
		observations = data[noNA,observations]
		
		
	} else {
		# observations must be a vector and covariates a matrix
		trend = as.matrix(trend)
		theNA = is.na(data) | apply(trend, 1, 
				function(qq) any(is.na(qq))
		)
		noNA = !theNA
		
		observations=data[noNA]
		covariates=trend[noNA,,drop=FALSE]
	}
	
	if(any(theNA)) {
		if(length(grep("^SpatialPoints", class(coordinates)))) {
			coordinates = SpatialPoints(coordinates)[noNA]	
		} else if(class(coordinates)=="dist"){
			coordinates = as.matrix(coordinates)
			coordinates = coordinates[noNA,noNA]
			coordinates = as.dist(coordinates)
			
		} else {
			warning("missing vlaues in data but unclear how to remove them from coordinates")
		}
	}
	
	# if the model's isotropic, calculate distance matrix
	
	
	if(!length(grep("^aniso", paramToEstimate)) &
			length(grep("^SpatialPoints", class(coordinates)))) {
		
		# see if there's anisotropy which is fixed, not estimated
		# which would be odd but we'll test nonetheless
		if(length(grep("^anisoRatio", names(param)))){
			if(abs(param["anosi.ratio"]- 1) > 0.0001){
				# it is indeed anisotropic
				
				if(any(names(param)=="anisoAngleDegrees") & 
						!any(names(param)=="anisoAngleRadians") ) {
					param["anisoAngleRadians"] = param["anisoAngleDegrees"]*2*pi/360				
				}
				
				x = coordinates@coords[,1] + 1i*coordinates@coords[,2]
				
				x = x * exp(-1i*param["anisoAngleRadians"])
				x = Re(x) +  (1i/ param["anisoRatio"] )*Im(x)
				coordinates = SpatialPoints(cbind(Re(x), Im(x)))
			} # end is anisotripic		
		} # end anisotropy params supplied
		coordinates = dist(coordinates@coords)		
		parscaleDefaults["range"] = sd(coordinates)/20
	} else {
		# doing geometric anisotropy
		parscaleDefaults["range"] = dist(t(bbox(coordinates)))/100
	}
	
	parscaleDefaults[names(parscale)] = parscale
	
	
	# default starting values for parameters
	paramDefaults = c(nugget=0,anisoRatio=0, anisoAngleDegrees=0,
			anisoAngleRadians=0,shape=1, boxcox=1,
			parscaleDefaults["range"])
	
	startingParam = param[paramToEstimate]
	names(startingParam) = paramToEstimate # fixes names lost when no starting value provided
	
	naStarting = is.na(startingParam)
	startingParam[naStarting]= paramDefaults[names(startingParam)[naStarting]]
	
	moreParams = param[!names(param) %in% paramToEstimate]
	
	
	# check to see if it's worth storing box cox quantities
	if(any(names(param)=="boxcox") & !any(paramToEstimate=="boxcox")) {
		if(abs(param["boxcox"]-1)>0.0001){ # boxcox not 1
			stored = list(
					boxcox = param["boxcox"],
					twoLogJacobian = 2*(param["boxcox"]-1)* 
							sum(log(observations))	
			)
			if(abs(param["boxcox"])<0.001) {
				stored$observations = log(observations) 
			} else  { #boxcox far from 0 and 1
				stored$observations <- 
						((observations^param["boxcox"]) - 1)/
						param["boxcox"]
			}
			
		} else { # boxcox is 1
			stored=NULL
		}
	} else { # no box cox
		stored=NULL
	}

	
#	observations, covariates, coordinates,





	
resG = (logLikGradient(param=startingParam,
					moreParams=moreParams, 
					observations=observations, 
					covariates=covariates, coordinates=coordinates,
					reml=reml, 
					eps=parscaleDefaults[paramToEstimate],
					lowerG=lowerDefaults[names(startingParam)], 
					upperG=upperDefaults[names(startingParam)]
			))
	
	fromOptim = optim(fn=loglikLgmG, par=startingParam, 
			gr=	logLikGradient,
			lower=lowerDefaults[names(startingParam)], 
			lowerG=lowerDefaults[names(startingParam)], 
			upper=upperDefaults[names(startingParam)],
			upperG=upperDefaults[names(startingParam)],
			eps=parscaleDefaults[names(startingParam)],
			observations=observations, 
			covariates=covariates, coordinates=coordinates,
			reml=reml, moreParams=moreParams,
			method = "L-BFGS-B", stored=stored
	)

return(fromOptim)
	
	fromLogLik = loglikLgm(param=fromOptim$par,
			moreParams=moreParams, 
			data=observations, 
			formula=covariates, coordinates=coordinates,
			reml=reml
	)
	
	result = list(
			param=c(attributes(fromLogLik)$betaHat,
					attributes(fromLogLik)$param),
			varBetaHat = attributes(fromLogLik)$varBetaHat,
			optim=fromOptim, 
			resid = as.vector(attributes(fromLogLik)$resid)
	)
	if(class(trend)=="formula") {
		result$trend = trend
	} else {
		result$formula= names(trend)
	}
	
	
	parameterTable = data.frame(estimate=result$param)
	rownames(parameterTable) =  names(result$param)
	
	parameterTable$stdErr = NA
	
	stdErr = sqrt(diag(attributes(fromLogLik)$varBetaHat))
	# sometimes varBetaHat doesn't have names
	parameterTable[names(attributes(fromLogLik)$betaHat), "stdErr"] =
			stdErr
	
	thelims = c(0.01, 0.99, 0.025, 0.975, 0.1, 0.9)
	for(D in thelims) {
		theQ = qchisq(D, 1, lower.tail=FALSE)
		parameterTable[[paste("ci",D,sep="")]] =
				parameterTable$estimate - sqrt(theQ) * parameterTable$stdErr
	}
	
	parameterTable[,"pr(est|par=0)"] = pchisq(
			parameterTable$estimate^2  / parameterTable$stdErr^2,
			df=1,lower.tail=FALSE)
	
	parameterTable[,"Estimated"] = FALSE
	parameterTable[paramToEstimate,"Estimated"] = TRUE
	parameterTable[names(attributes(fromLogLik)$betaHat),"Estimated"] = TRUE
	if(estimateVariance)
		parameterTable["variance","Estimated"] = TRUE
	
	
	rownames(parameterTable)=gsub("^variance$", "sdSpatial", 
			rownames(parameterTable))	
	rownames(parameterTable)=gsub("^nugget$", "sdNugget", 
			rownames(parameterTable))	
	parameterTable[c("sdSpatial", "sdNugget"),"estimate"] = 
			sqrt(parameterTable[c("sdSpatial", "sdNugget"),"estimate"])
	
	result$summary = parameterTable
	
	result
}


