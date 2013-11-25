loglikLgm = function(param,
		data, trend, coordinates=data,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=NULL) {

	# create 'covariates', and 'observations'
 	
	if(class(trend)=="formula") {
		covariates = model.matrix(trend, as.data.frame(data))
		observations = formulaLhs(trend)
		if(!any(names(data)==observations))
			warning("can't find observations ", observations, "in data")
		observations = data[[observations]]
	} else {
		# observations must be a vector and covariates a matrix
		observations=data
		covariates=as.matrix(trend)
	}
	
		if(!length(grep("SpatialPoints", class(coordinates))) &
				class(coordinates) != "dist")	{
			warning("coordinates must be a SpatialPoints object\n or a dist object.  It's being assumed \n it's a matrix of coordinates")
			coordinates = SpatialPoints(coordinates)
		}

		# sort out the parameters
		param = c(param, moreParams)
		
		names(param) = gsub("^var$", "variance", names(param))
		
		param=param[!is.na(param)]
	
		haveVariance = any(names(param)=="variance")
		haveNugget = any(names(param)=="nugget")
		if(haveNugget){
			if(param["nugget"] == 0) {
				haveNugget=FALSE
			}
		} else {
			param["nugget"] = 0
		}

		if(!haveVariance) {
			# variance will be estimated.  
		  # nugget is assumed to be nugget/variance
			param["variance"] = 1
		}
			
		# box cox transform
		twoLogJacobian = 0
		if(any(names(param)=="boxcox")) {
			
			if(any(names(stored)=="boxcox")) {
				if(abs(stored$boxcox-param["boxcox"])>0.0001)
					warning("boxcox param and stored are ", param["boxcox"],
							" and ", stored$boxcox)
				twoLogJacobian = stored$twoLogJacobian
				observations = stored$observations
			} else { # need to compute boxcox stuff
			
			if(abs(param["boxcox"]-  1 ) < 0.001) {
			 # boxcox close to 1, don't transform
				twoLogJacobian=0
				
			} else { # box cox is not one.
				twoLogJacobian = 2*(param["boxcox"]-1)* 
					sum(log(observations))	
			
				if(is.nan(twoLogJacobian))
					warning("boxcox shouldnt be used with negative data")

				if(abs(param["boxcox"])<0.001) {
					observations = log(observations) 
				} else  { #boxcox far from 0 and 1
				observations <- 
						((observations^param["boxcox"]) - 1)/
							param["boxcox"]
				}
			} # end boxcox param far from 1
		} # end need to compute boxcox
		} # end have box cox
		
		
		# variance matrix
		covMat = matern(x=coordinates, param=param)	

		if(haveNugget)
			Matrix::diag(covMat) = Matrix::diag(covMat) + param["nugget"]

		# matrix operations
		cholCovMat = Matrix::chol(covMat)
	
		# cholCovMat %*% t(cholCovMat) = covMat
		cholCovInvX = Matrix::solve(cholCovMat, covariates)
		cholCovInvXcross = Matrix::crossprod(cholCovInvX)
	
		cholCovInvXcrossInv = Matrix::solve(cholCovInvXcross)
	
		# cholCovInvX = cholCovMat^{-1} %*% covariates
		cholCovInvY = Matrix::solve(cholCovMat, observations)
		# cholCovInvY = cholCovMat^{-1} %*% observations
		
		# beta hat = (D'V^{-1}D)^{-1}D'V^{-1}y
		# V = L L', V^{-1} = t(L^(-1)) %*% L^(-1)
		# beta hat = (D' Linv' Linv D)^{-1}D'Linv' L y

		# Covariates and likelihood
	
		betaHat = cholCovInvXcrossInv %*% Matrix::crossprod(cholCovInvX, cholCovInvY) 
		
		resids = observations - covariates %*% betaHat
		# sigsqhat = resids' %*% Vinv %*% residsx
		#    =   resids' Linv' Linv resids
		cholCovInvResid = Matrix::solve(cholCovMat, resids)
		
		totalSsq = as.vector(Matrix::crossprod(cholCovInvResid))

		totalVarHat = totalSsq/length(observations)
		
		if(!haveVariance) { # profile likelihood with optimal sigma
		
			minusTwoLogLik = length(observations) * log(2*pi) + 
				length(observations) * log(totalVarHat) +
				2*sum(log(Matrix::diag(cholCovMat))) +
				length(observations) - twoLogJacobian
		
			if( reml ) {
			# REML
				choldet =  2*sum(log(Matrix::diag(chol(cholCovInvXcross)))) 
		
				minusTwoLogLik =  minusTwoLogLik + choldet
		
				}  # end reml
			
		param[c("variance","nugget")] = totalVarHat * param[c("variance","nugget")]

	} else { # a variance was supplied
		# calculate likelihood with the variance supplied
		# -2 log lik = n log(2pi) + log(|V|) + resid' Vinv resid
		
		minusTwoLogLik = length(observations) * log(2*pi) +
				2*sum(log(Matrix::diag(cholCovMat))) +
					totalSsq - twoLogJacobian
		totalVarHat = 1
		if( reml ) {
			reml=FALSE
			warning("a variance was supplied by REML was requested.  doing ML instead.")		
		}
	}

	# format the output
	betaHat = as.vector(betaHat)
	names(betaHat) = colnames(covariates)
	varBetaHat = totalVarHat * cholCovInvXcrossInv 

	result = minusTwoLogLik
	if(minustwotimes) {
			names(result) = "minusTwoLogLik"
	} else {
			result = -0.5*result

			names(result)="logLik"
	} 
	if(reml)
		names(result) = gsub("Lik$", "RestrictedLik", names(result))
	
	attributes(result)$param = param
#		attributes(result)$totalVarHat = totalVarHat
	attributes(result)$betaHat = betaHat
	attributes(result)$varBetaHat = as.matrix(varBetaHat)
#	attributes(result)$reml=reml
#		attributes(result)$twoLogJacobian = twoLogJacobian
#		attributes(result)$choldet = as.vector(choldet)
	attributes(result)$resid = resids
	result
}

 


likfitLgm = function(
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
	lowerDefaults = c(nugget=0,range=0,anisoRatioratio=0.0001,
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
	
	fromOptim = optim(fn=loglikLgm, par=startingParam, 
			lower=lowerDefaults[paramToEstimate], 
			upper=upperDefaults[paramToEstimate],
			control = list(parscale=parscaleDefaults[paramToEstimate]),
			data=observations, 
			trend=covariates, coordinates=coordinates,
		reml=reml, moreParams=moreParams,
		method = "L-BFGS-B", stored=stored
		)
	
		
	fromLogLik = loglikLgm(param=fromOptim$par,
			moreParams=moreParams, 
			data=observations, 
			trend=covariates, coordinates=coordinates,
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
		result$trend= names(trend)
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
	
	dimnames(parameterTable) = lapply(dimnames(parameterTable),
			function(qq) {
				qq=gsub("_", "\\\\textunderscore ", qq)
				qq=gsub("\\$", "\\\\textdollar ", qq)
				qq=gsub("<", "\\\\textless ", qq)
				qq=gsub(">", "\\\\textgreater ", qq)
				qq
			}
	)
	
	
	result$summary = as.data.frame(parameterTable)
	
	result
}


