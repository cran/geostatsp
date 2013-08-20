
loglikP = function(param, 
		observations, covariates, coordinates,
		method=c("ML","REML"), paramNames = names(param)) {
	
		names(param) = paramNames
		names(param) = gsub("^var$", "variance", names(param))
		
	
		haveVariance = any(names(param)=="variance")
		haveNugget = any(names(param)=="nugget")
		
		if(!haveVariance) {
			# variance will be estimated.  
		  # nugget is assumed to be nugget/variance
			param["variance"] = 1
		}
			
		if(!haveNugget)
			param["nugget"] = 0
		
		covMat = matern(x=coordinates, param=param)		
		
		
		cholCovMat = chol(covMat)
		# cholCovMat %*% t(cholCovMat) = covMat
		cholCovInvX = solve(cholCovMat, covariates)
		# cholCovInvX = cholCovMat^{-1} %*% covariates
		cholCovInvY = solve(cholCovMat, observations)
		# cholCovInvY = cholCovMat^{-1} %*% observations
		
		# beta hat = (D'V^{-1}D)^{-1}D'V^{-1}y
		# V = L L', V^{-1} = t(L^(-1)) %*% L^(-1)
		# beta hat = (D' Linv' Linv D)^{-1}D'Linv' L y
	
		
 		
}


if(F) {
	
	observations = rnorm(10)
	covariates = cbind(rep(1, length(observations)), 1:length(observations))
	coordinates =  SpatialPointsDataFrame(
			cbind(runif(length(observations)), 
					runif(length(observations))),
			data=data.frame(id=1:length(observations)))
	method="ML"
	
	param = c(range=2, rough=2,	nugget=0.3,
			aniso.ratio=0.6, aniso.angle.degrees=30)
	paramNames = names(param)
	
	
}