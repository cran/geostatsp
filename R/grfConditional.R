grfConditional = function(data,  model.fit, locations, Nsim, covariates,rasterMethod = "ngb", 
		fun, nugget.in.prediction=T ) {
	
	# mean of obseved data
	xmat = model.matrix(model.fit$formula, data@data)
	if(dim(xmat)[1] < dim(data@data)[1]) {
		# some missing values
		notNA = rownames(data@data) %in% rownames(xmat)
		data = data[notNA,]
		xmat = model.matrix(model.fit$formula, data@data)
	}
	
	thefit = xmat %*% model.fit$beta[colnames(xmat)]
	theterms = strsplit(as.character(model.fit$formula), "~")
	U = data@data[,theterms[[2]]] -thefit
	
	
	# find grid on which predictions are made for fixed effects	
	if(!missing(covariates)) {
	if(is.list(covariates)){
		if(length(rasterMethod)!= length(covariates)) {
			themethod = rep(rasterMethod[1], length(covariates))
			names(themethod) = names(covariates)
		}
		locations.mean = stackRasterList(covariates,locations, themethod)	
	}  else {
		locations.mean = covariates
	}
	} # end if missing covariates
	if(missing(covariates) & !missing(locations))
		locations.mean = locations
	
	
	if(missing(locations) & !missing(covariates))
		locations = raster(covariates)

	if(!missing(locations) & !missing(covariates)){
		if(!is.null(covariates))
			locations.mean = stackRasterList(list(locations.mean),template=locations)
	} 
	
	if(missing(covariates) & missing(locations))
		warning("one of covariates and locatoins must be supplied")
	
	
	coordsMat = as.data.frame(locations, xy=T)[,c("x","y")]
	colnames(data@coords)=c("x","y")
	allCoords = as.matrix(rbind(coordsMat, data@coords))
 	
	anisoMat =
			matrix(c(cos(model.fit$aniso.pars["psiA"]), sin(model.fit$aniso.pars["psiA"]), 
							-sin(model.fit$aniso.pars["psiA"]), cos(model.fit$aniso.pars["psiA"])),2) %*% 
			diag(c(1/(model.fit$phi*model.fit$aniso.pars["psiR"]),1/model.fit$phi))
	
	
	coordsTrans = allCoords %*% anisoMat
	
	distmat = as.matrix(dist(coordsTrans))
	corMat = model.fit$sigmasq*geoR::matern(distmat,1,model.fit$kappa)
	
	
	Npred = dim(coordsMat)[1]
	Spred = 1:Npred
	Ndata=length(data)
	Ncoords = dim(allCoords)[1]
	Sdata = seq(Npred+1, Ncoords)
	sig11 = corMat[Spred,Spred]
	if(nugget.in.prediction)
		diag(sig11) = diag(sig11) + model.fit$nugget
	sig12 = corMat[Spred, Sdata]
	sig22 = corMat[Sdata, Sdata]
	diag(sig22) = diag(sig22)+model.fit$nugget
	sig22Inv = chol2inv(chol(sig22))
	
	condVar = sig11 + sig12 %*% sig22Inv %*% t(sig12)
	condVarChol = as.matrix(chol(condVar))

	
	
	
	
	
	pred = sig12 %*% sig22Inv %*% U
	
	
	# grid on which spatial predictions are made	 
	
	
	if(!missing(covariates)) {
		if(dim(xmat)[2]==2)
			names(locations.mean) = colnames(xmat)[2]
		covariates = as.data.frame(locations.mean)
		covariates[,theterms[[2]]]=0
		theNA = apply(covariates, 1, function(qq) any(is.na(qq)))
		covariates[is.na(covariates)] = 0
		xmat = 	 model.matrix(model.fit$formula, covariates)
		thefitPred = xmat %*% model.fit$beta[colnames(xmat)]
		thefitPred[theNA] = NA
		predPLusFit = pred +thefitPred
	} else {
		predPLusFit = pred
	}	
	
	

	
	if(missing(fun)) {
		simBig=rep(0.1, Nsim*Npred)
		for(D in 1:Nsim) {
			sim = predPLusFit + condVarChol %*% rnorm(Npred)
 
			simBig[seq(Npred*(D-1)+1, len=Npred)] = sim
		}

	result = brick(array(0.1,c(locations@nrows, locations@ncols,Nsim)), xmn=locations@extent@xmin, 
			xmx=locations@extent@xmax, ymn=locations@extent@ymin, 
			ymx=locations@extent@ymax,crs=locations@crs )
	values(result) = simBig	
	} else { # a function was supplied

		
result = list()			
		for(D in 1:Nsim) {
			result[[D]] = fun(predPLusFit + condVarChol %*% rnorm(Npred))			
		}
	}
	result
}

