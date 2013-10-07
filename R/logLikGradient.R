
logLikGradient = function(param, 
		observations, covariates, coordinates,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=NULL,
		upper, lower, eps=rep(0.001, length(param))
	){
		
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
	


	# the function to get gradients
grFun = function(parname, right=TRUE) {
	
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

	# calcluate the gradient
	parForMapply = rep(names(param), 2)
	rightForMapply = rep(c(TRUE, FALSE), rep(length(param),2))

	thegrad = mapply(grFun, parname=parForMapply, right=rightForMapply)
	
	thegrad = matrix(thegrad, nrow=length(names(param)), ncol=2, 
					dimnames = list(names(param), 
							c("left","right")[1+rightForMapply[c(1,length(param)+1)]]))
			
	thegrad = thegrad[,"right"] - thegrad[,"left"]
	thegrad = thegrad * (paramsRight[names(thegrad)] - paramsLeft[names(thegrad)])
	thegrad
	
}

if(FALSE) {
	
	
	n=100
	mydat = SpatialPointsDataFrame(cbind(runif(n), runif(n)), 
			data=data.frame(cov1 = rnorm(n), cov2 = rpois(n, 0.5))
	)
	
	trueParamAniso = param=c(variance=2^2, range=0.2, rough=2,
			nugget=0,aniso.ratio=4,aniso.angle.degrees=10, nugget=0)
	
	mydat$U = GaussRF(mydat, par=trueParamAniso)
	mydat$Y = 12 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
			mydat$U + rnorm(length(mydat), 0, sd=sqrt(trueParamAniso["nugget"]))
	
	
	myres = likfitLgm(mydat, Y ~ cov1 + cov2, 
			param=c(range=0.1,nugget=0.1,rough=2), 
			paramToEstimate = c("range","nugget", "aniso.ratio","aniso.angle.degrees")
	)
	
	trendmat = cbind(1,as.matrix(mydat@data[,c("cov1","cov2")]))
	
	temp=loglikLgm(param=myres$param, 
			data=mydat$Y, trend=trendmat, coordinates=mydat,
			reml=FALSE, minustwotimes=TRUE)
	
	
	param = myres$param[c("range","nugget","aniso.ratio","aniso.angle.degrees","rough")]
	lower = c(range=0,nugget=0,aniso.ratio = 0, aniso.angle.degrees=-90, rough=0)
	upper = c(range=Inf,nugget=Inf,aniso.ratio = Inf, aniso.angle.degrees=90, rough=Inf)
	
			logLikGradient( param, 
				mydat$Y, trendmat, mydat,
					reml=TRUE, 
					minustwotimes=TRUE,
					stored=NULL, moreParams=NULL,
					upper=upper, lower=lower, eps=rep(0.001, length( param))
		)
	
logLikGradient( param, 
		mydat$Y, trendmat, mydat,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=c(boxcox=0.5),
		upper=upper, lower=lower, eps=rep(0.001, length( param))
)
logLikGradient( param, 
		mydat$Y, trendmat, mydat,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=c(boxcox=0),
		upper=upper, lower=lower, eps=rep(0.001, length( param))
)

param2 = c(param, boxcox=1)
lower2 = c(lower, boxcox=-Inf)
upper2 = c(upper, boxcox=Inf)
logLikGradient( param2, 
		mydat$Y, trendmat, mydat,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=NULL,
		upper=upper2, lower=lower2, eps=rep(0.001, length( param2))
)

toFix = c("aniso.angle.degrees","rough")
param2 = param[!names(param) %in% toFix]
upper2 = upper[!names(upper) %in% toFix]
lower2 = lower[!names(lower) %in% toFix]
logLikGradient( param2, 
		mydat$Y, trendmat, mydat,
		reml=TRUE, 
		minustwotimes=TRUE,
		stored=NULL, moreParams=param[toFix],
		upper=upper2, lower=lower2, eps=rep(0.001, length( param2))
)

}
