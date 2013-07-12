
loglik.GRF = function(geodata, ...) {
	UseMethod("loglik.GRF")
	
}


loglik.GRF.SpatialPointsDataFrame <- function(geodata, 
		formula, sigmasq, phi,
	psiA=1,psiR=1,kappa=1,lambda=1,
	nugget=0,  ...){

 
	theterms = strsplit(as.character(formula), "~")
	
	theCovs = attributes(terms(formula))$term.labels
	
	newdata = as.geodata(geodata, 
			data.col=theterms[[2]], 
			covar.col=theCovs)

	trend = as.formula(paste("~", theterms[[3]]))

	# for some reason this line is needed!
	cov.model = list(...)$cov.model

	if( !any( cov.model == 
		c("matern", "exponential", "gaussian", "spherical", "circular", "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", "gneiting", "gneiting.matern", "pure.nugget")
		)
	) warning("cov.model not implemented in geoR")


result=loglik.GRF(geodata = newdata,
		cov.model=cov.model,
		cov.pars=c(sigmasq, phi),
		nugget=nugget, kappa=kappa,
		lambda=lambda, psiR= psiR,
		psiA=psiA, trend= trend,
		method.lik="ML",
		compute.dists=TRUE)


 

result
}



loglik.GRF.default = 
		function(geodata, coords=geodata$coords, data=geodata$data,
				obj.model = NULL,
				cov.model="exp", cov.pars,
				nugget=0, kappa=0.5, lambda=1, psiR=1, psiA=0,
				trend="cte", method.lik="ML",
				compute.dists = TRUE, realisations = NULL, ...)
{
	if(!is.null(obj.model)){
		if(!is.null(obj.model$cov.model)) cov.model <- obj.model$cov.model
		if(!is.null(obj.model$cov.pars)) cov.pars <- obj.model$cov.pars
		if(!is.null(obj.model$nugget)) nugget <- obj.model$nugget
		if(!is.null(obj.model$kappa)) kappa <- obj.model$kappa
		if(!is.null(obj.model$lambda)) lambda <- obj.model$lambda
		if(!is.null(obj.model$psiR)) psiR <- obj.model$psiR
		if(!is.null(obj.model$psiA)) psiA <- obj.model$psiA      
		if(!is.null(obj.model$trend)) trend <- eval(obj.model$trend)
		## a resolver: problema em passando  trend
	}
	sigmasq <- cov.pars[1]
	phi <- cov.pars[2]
	if(method.lik == "REML" | method.lik == "reml" | method.lik == "rml") 
		method.lik <- "RML"
	if(method.lik == "ML" | method.lik == "ml")
		method.lik <- "ML"
	if(is.null(realisations))
		realisations <- as.factor(rep(1, length(data)))
	else
		realisations <- as.factor(realisations)
	nrep <- length(levels(realisations))
	##
	## Absurd values
	##
	if(kappa < 1e-04) return(-(.Machine$double.xmax^0.5))
	if((nugget+sigmasq) < 1e-16) return(-(.Machine$double.xmax^0.5))
	##
	## Trend matrix
	##
	if(missing(geodata))
		xmat <- trend.spatial(trend=trend, geodata=list(coords = coords, data = data))
	else xmat <- unclass(trend.spatial(trend=trend, geodata=geodata))
	if (nrow(xmat) != nrow(coords)) 
		stop("coords and trend have incompatible sizes")
	beta.size <- ncol(xmat)
	xmat <- split(as.data.frame(unclass(xmat)), realisations)
	xmat <- lapply(xmat, as.matrix)
	##
	## Anisotropy
	##
	vecdist <- function(x){as.vector(dist(x))}
	if(psiR != 1 | psiA != 0){
		coords.c <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
		.likGRF.dists.vec <- lapply(split(as.data.frame(coords.c),
						as.factor(realisations)), vecdist)
	}
	else if(compute.dists) .likGRF.dists.vec <- lapply(split(as.data.frame(coords),
						as.factor(realisations)), vecdist)
	##
	## Box-Cox transformation
	##
	z <- data
	if(abs(lambda - 1) < 0.0001)
		log.jacobian <- 0
	else {
		if(any(z <= 0))
			stop("Transformation not allowed for zero or negative data")
		data <- z^(lambda - 1)
		if(any(data <= 0)) log.jacobian <- log(prod(data))
		else log.jacobian <- sum(log(data))
		data <- NULL
		if(abs(lambda) < 0.0001)
			data <- log(z)
		else data <- ((z^lambda) - 1)/lambda
	}
	data <- split(data, as.factor(realisations))
	##
	## Computing likelihood
	##
	sumnegloglik <- 0
	for(i in 1:nrep){
		## NOTE: Likelihood for Independent observations 
		##       arbitrary criteria used here:
		##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
		##
		n <- length(data[[1]])
		if((phi < 1e-16) | (sigmasq < 1e-16)){
			V <- list(varcov = diag(x=(nugget+sigmasq), n),
					log.det.to.half = (n/2) * log(nugget+sigmasq))
		}
		else{
			V <- varcov.spatial(dists.lowertri = .likGRF.dists.vec[[i]],
					cov.model = cov.model, kappa=kappa,
					nugget = nugget, cov.pars=c(sigmasq, phi),
					det = TRUE)
		}
		if(!is.null(V$crash.parms)){
			cat("varcov.spatial: improper matrix for following the given parameters:")
			print(V$crash.parms)
			stop()
		}
		ivx <- solve(V$varcov,xmat[[i]])
		xivx <- crossprod(ivx,xmat[[i]])
		betahat <- geoR:::.solve.geoR(xivx, crossprod(ivx,data[[i]]))
		res <- data[[i]] - xmat[[i]] %*% betahat
		ssres <- drop(crossprod(res, solve(V$varcov,res)))
		if(method.lik == "ML"){
			negloglik <- (n/2)*(log(2*pi)) + V$log.det.to.half +  0.5 * ssres
		}
		if(method.lik == "RML"){
			choldet <- sum(log(diag(chol(xivx))))
			negloglik <- V$log.det.to.half +  0.5 * ssres + choldet
			xx.eigen <- eigen(crossprod(xmat[[i]]), symmetric = TRUE, only.values = TRUE)
			negloglik <- negloglik + ((n-beta.size)/2)*(log(2*pi)) - 0.5 * sum(log(xx.eigen$values))
		}
		sumnegloglik <- sumnegloglik + negloglik
	}
	sumnegloglik <- sumnegloglik - log.jacobian
	if(sumnegloglik > (.Machine$double.xmax^0.5))
		sumnegloglik <- .Machine$double.xmax^0.5
	return(as.vector(-sumnegloglik))
}
