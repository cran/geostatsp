profLlgm = function(fit,mc.cores=NULL, ...) {
	
	dots = list(...)
	varying = intersect(names(dots), names(fit$param))

	nonLinearParams = c('boxcox','shape','nugget','variance',
			'anisoAngleDegrees','anisoRatio','range')
	
	reEstimate = rownames(fit$summary)[
			fit$summary[,"Estimated"]
	]
	reEstimate = gsub("sdNugget", "nugget", reEstimate)
	reEstimate = gsub("sdSpatial", "variance", reEstimate)
	reEstimate = gsub("range \\(km\\)", "range", reEstimate)
	reEstimate = intersect(reEstimate, nonLinearParams)
	reEstimate = reEstimate[!reEstimate %in% varying]
	
	parValues = do.call(expand.grid, dots[varying])
	
	baseParams = fit$param
	baseParams = baseParams[names(baseParams)%in%
					nonLinearParams]
	
	baseParams=baseParams[!names(baseParams)%in% varying]
	baseParams=baseParams[names(baseParams) != 'variance']
	
	
	
	oneL = function(...){
		res=geostatsp::likfitLgm(fit$data, 
				fit$model$trend,
				param=c(..., 
						baseParams),
				paramToEstimate=reEstimate,
				reml=fit$model$reml)
		c(..., L = 	res$opt$val,
				res$param[reEstimate])
	}
	
	forCall = c(
			as.list(parValues),
			FUN=oneL,
			SIMPLIFY=TRUE)
	
	if(!is.null(mc.cores)) {
		resL = do.call(parallel::mcmapply,
				c(forCall, mc.cores=mc.cores)
		)
	} else {
		resL = do.call(mapply,forCall)
	}

	
	if(length(varying)==1) {
		dots[[varying]] = c(
				dots[[varying]],
				fit$param[varying]
		)
		theorder = order(dots[[varying]])
		L = c(resL["L",],fit$opt$val)[theorder]
		dots[[varying]] = dots[[varying]][theorder]
		
	} else {
		thedimnames=dots[varying]
		for(D in names(thedimnames))
			thedimnames[[D]] = paste(
					D, thedimnames[[D]],sep='_')
		L = array(resL["L",], 
				unlist(lapply(thedimnames, length)),
				dimnames=thedimnames)	
	} 
	
	Sprob = c(1, 0.995, 0.99, 0.95, 0.9, 0.8, 0.5, 0)
	Squant = qchisq(Sprob, df=length(varying))
	Scontour = -fit$opt$val/2 -Squant/2 	
	
	
	res = list(logL=-L/2,
			full=t(resL),
			prob=Sprob,
			breaks=Scontour,
			MLE=fit$param[varying],
			maxLogL = -fit$opt$val/2,
			basepars=baseParams
	)


#	dput( rev(
#			RColorBrewer::brewer.pal(
#					length(Scontour)-1,
#					'Spectral')	))
	res$col=c("#3288BD", "#99D594", 
			"#E6F598", "#FFFFBF", "#FEE08B", "#FC8D59", 
			"#D53E4F")
	
	
	res$breaks[1] = res$breaks[2]-abs(min(res$logL))
	names(res$col) = as.character(res$prob[-length(res$prob)])
	
	res = c(dots[varying],res)
	
	if(length(varying)==1) {
		smaller = dots[[varying]] <= res$MLE
		bigger = dots[[varying]] >= res$MLE
		Skeep = seq(2, length(Sprob)-1)
		res$ci = cbind(
				prob=Sprob[Skeep],
				lower= approx(res$logL[smaller], 
						dots[[varying]][smaller],
						Scontour[Skeep])$y,
				upper=approx(res$logL[bigger], 
						dots[[varying]][bigger], 
						Scontour[Skeep])$y
		)
		
		res$ciLong = na.omit(
				reshape(as.data.frame(res$ci), 
						direction="long",
						varying=list(par=c('upper','lower')),
						v.names='par',
						times = c('upper','lower'),
						timevar=c('direction'),
						idvar='prob')	
		)
		res$ciLong = rbind(
				res$ciLong,
				data.frame(prob=0,direction='upper',
						par=res$MLE)[
						colnames(res$ciLong)
				])
		res$ciLong$quantile = (1-res$ciLong$prob)/2
		res$ciLong[res$ciLong$direction=='upper','quantile'] =
				1 - res$ciLong[res$ciLong$direction=='upper','quantile'] 
		
		
		res$ciLong = res$ciLong[order(res$ciLong$par),]
	}
	
	res
}
