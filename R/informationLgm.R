

informationLgm = function(fit, ...) {
	nonLinearParams = c('boxcox','shape','nugget','variance',
			'anisoAngleDegrees','anisoRatio','range')
	
	reEstimate = rownames(fit$summary)[
			fit$summary[,"Estimated"]
	]
	reEstimate = gsub("sdNugget", "nugget", reEstimate)
	reEstimate = gsub("sdSpatial", "variance", reEstimate)
	reEstimate = gsub("range \\(km\\)", "range", reEstimate)
	reEstimate = intersect(reEstimate, nonLinearParams)
	
	baseParam = fit$param[reEstimate]
	moreParams = fit$param[
			!names(fit$param) %in% reEstimate &
					names(fit$param) %in% nonLinearParams]
	
	parToLog = c("nugget","variance","anisoRatio","range")
	parToLog = intersect(reEstimate, parToLog)
	
	if(!all(baseParam[parToLog]>0))
		return(list(summary=fit$summary,information=NULL))

	
	
	
	oneL = function(param, ...) {
		param[parToLog] = exp(param[parToLog])
		loglikLgm(param, ...)
	}
	
	baseParam[parToLog] = log(baseParam[parToLog])
	
	# get rid of NA's
	fit$data = fit$data[,]
	
	hess = numDeriv::hessian(oneL, baseParam,
			data=fit$data,formula=fit$model$formula,
			reml=fit$model$reml,
			moreParams=moreParams, ...)
	
	whichLogged = which(names(baseParam)%in% parToLog)
	names(baseParam)[whichLogged] = paste("log(", 
			names(baseParam)[whichLogged], ")",sep="")
	
	dimnames(hess) = list(names(baseParam),names(baseParam))
	infmat = solve(hess)*2
	
	pvec = grep("^ci([[:digit:]]|\\.)+$", colnames(fit$summary),
			value=TRUE)
	pvec = as.numeric(gsub("^ci","", pvec))
	qvec = qnorm(pvec)
	names(qvec) = paste("ci", pvec, sep="")
	
	stdErr = diag(infmat)
	if(any(stdErr<0))
		return(list(summary=fit$summary,information=infmat))
	stdErr = sqrt(stdErr)
	
	toAdd = outer(stdErr, qvec, FUN="*")

	forSummary = baseParam[rownames(toAdd)] + toAdd
	expAdd = exp(forSummary[
					grep("^log\\(",rownames(forSummary))
					,])
	rownames(expAdd) = gsub("^log\\(|\\)$","",rownames(expAdd))
	forSummary = rbind(forSummary, expAdd)	

	summary = fit$summary
	
	if(any(rownames(forSummary)=='nugget'))
		forSummary = rbind(forSummary,
				sdNugget = sqrt(
						pmax(0,forSummary["nugget",]))
		)
	if(any(rownames(forSummary)=='variance'))
		forSummary = rbind(forSummary,
				sdSpatial = sqrt(
						pmax(0,forSummary["variance",]))
		)
	
	
	inBoth = intersect(rownames(summary), rownames(forSummary))
	
	summary[inBoth,colnames(forSummary)] = 
			forSummary[inBoth,]
	
	list(summary=summary,information=infmat)
}
