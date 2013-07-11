likfit = function(geodata, ...) {
	UseMethod("likfit")
	
}
likfit.SpatialPointsDataFrame <- function(geodata, 
		formula, dist.rel.scale=100,  ...){


	
	theterms = unlist(strsplit(as.character(formula), "~"))
	theterms = gsub("[[:space::]", "", theterms)
	theterms = theterms[theterms != ""]
	
	theCovs = attributes(terms(formula))$term.labels
	
	newdata = as.geodata(geodata, 
			data.col=theterms[1], 
			covar.col=theCovs)

	trend = as.formula(paste("~", theterms[2]))
	
	
	
	maxdist = dist(t(geodata@bbox))
	if(!"control" %in% names(list(...))) {
		control=list(parscale=c(maxdist/dist.rel.scale,0.1))
	}

	# for some reason this line is needed!
	cov.model = list(...)$cov.model

	if( !any( cov.model == 
		c("matern", "exponential", "gaussian", "spherical", "circular", "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", "gneiting", "gneiting.matern", "pure.nugget")
		)
	) warning("cov.model not implemented in geoR")

	result = likfit.default(geodata=newdata, trend=trend, control=control, ...)
	result$formula = formula
	
	
	result
}

"likfit.default" <-
		function (geodata, coords=geodata$coords, data=geodata$data,
				trend = ~1, ini.cov.pars,
				fix.nugget = FALSE, nugget = 0, 
				fix.kappa = TRUE, kappa = 0.5, 
				fix.lambda = TRUE, lambda = 1, 
				fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, 
				cov.model, realisations, lik.method = "ML",
				components = TRUE, nospatial = TRUE, limits = pars.limits(), 
				print.pars = FALSE, messages, ...) 
{

	name.geodata <- deparse(substitute(geodata))
	##
	## Checking input
	##
	call.fc <- match.call()
	ldots <- list(...)
	.likGRF.dists.vec=NULL

	
	temp.list <- list()
	temp.list$print.pars <- print.pars
	if(missing(messages))
		messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
	else messages.screen <- messages
	##
	if(!missing(ini.cov.pars)){
		if(any(class(ini.cov.pars) == "eyefit")){
			ini.cov.pars <- ini.cov.pars[[1]]
#      cov.model <- ini.cov.pars[[1]]$cov.model
#      kappa <- ini.cov.pars[[1]]$kappa
		}
		if(any(class(ini.cov.pars) == "variomodel")){
			cov.model <- ini.cov.pars$cov.model
			kappa <- ini.cov.pars$kappa
		}
	} else {
		# default initial values
		# 
		ini.cov.pars = c(var(data), sd(coords[,1]))	
	}

	
	if(missing(cov.model)) cov.model <- "matern"
	cov.model <- match.arg(cov.model, choices = geoR:::.geoR.cov.models)
	if(cov.model == "stable") cov.model <- "powered.exponential"
	if(any(cov.model == c("power", "gneiting.matern", "gencauchy")))
		stop(paste("parameter estimation for", cov.model, "is not yet implemented"))
	##  if(any(cov.model == c("gneiting.matern", "gencauchy"))){
	##    if(length(kappa != 2))
	##      stop(paste(cov.model, "requires two values in the argument kappa"))
	##    if(length(fix.kappa) == 1) fix.kappa <- rep(fix.kappa, 2) 
	##    stop("parameter estimation for gneiting.matern model is not yet implemented")
	##  }
	fixed.pars <- list(cov.model=cov.model)
	if(fix.nugget) fixed.pars$nugget <- nugget
	if(fix.kappa) fixed.pars$kappa <- kappa
	if(fix.psiA) fixed.pars$psiA <- psiA
	if(fix.psiR) fixed.pars$psiR <- psiR
	geoR:::.check.geoRparameters.values(list=fixed.pars, messages = messages.screen)
	if(cov.model == "matern" & all(kappa == 0.5)) cov.model <- "exponential"
	temp.list$cov.model <- cov.model
	if(cov.model == "powered.exponential")
		if(limits$kappa["upper"] > 2) limits$kappa["upper"] <- 2
	if(cov.model == "gencauchy")
		if(limits$kappa2["upper"] > 2) limits$kappa2["upper"] <- 2
	##
	## Likelihood method
	##
	#####
	##### temporary code back compatibility for argument "method"
	lik.MET <- c("ML", "ml", "RML", "REML", "rml", "reml")
	MET <- pmatch(names(ldots), "method") == 1
	if(!is.na(MET) && any(MET) && (ldots[[which(MET)]] %in% lik.MET)){
		warning("argument \"method\" has changed and is now used as an argument to be passed to optim(). Use \"lik.method\" to define the likelihood method")
		lik.method <- lik.MET[pmatch(ldots[[which(MET)]], lik.MET)]
		ldots[which(as.logical(pmatch(names(ldots), "method", nomatch=0)))] <- NULL
	}
	#####
	method.lik <- lik.method
	if(method.lik %in% c("REML","reml","rml","RML"))  method.lik <- "RML"
	if(method.lik %in% c("ML", "ml")) method.lik <- "ML"
	if(method.lik == "ML" & cov.model == "power")
		stop("\n\"power\" model can only be used with method.lik=\"RML\".\nBe sure that what you want is not \"powered.exponential\"")
	temp.list$method.lik <- method.lik
	##
	## setting coordinates, data and covariate matrices
	##

	coords <- as.matrix(coords)
	data <- as.vector(data)
	n <- length(data)

	
	if((nrow(coords) != n) | (2*n) != length(coords))
		stop("\nnumber of locations does not match with number of data")
	if(missing(geodata))
		xmat <- trend.spatial(trend=trend, geodata=list(coords = coords, data = data))
	else xmat <- unclass(trend.spatial(trend=trend, geodata=geodata))
	xmat.contrasts  <- attr(xmat,"contrasts")
	xmat <- unclass(xmat)
	if(nrow(xmat) != n)
		stop("trend matrix has dimension incompatible with the data")
	geoR:::.solve.geoR(crossprod(xmat))
	beta.size <- temp.list$beta.size <- dim(xmat)[2]
	##
	## setting a factor to indicate different realisations
	##
	
	if(missing(realisations))
		realisations <- as.factor(rep(1, n))
	else{
		if(!missing(geodata)){
			real.name <- deparse(substitute(realisations))
			if(all(isTRUE(as.logical(real.name))))
				if(is.null(geodata$realisations)) stop("element realisation not available in the geodata object")
				else realisations <- geodata$realisations
			else{
				if(!is.null(geodata[[real.name]]))
					realisations <- geodata[[real.name]]
			}
		}
		if(length(realisations) != n)
			stop("realisations must be a vector with the same length of the data")
		realisations <- as.factor(realisations)
	}
	temp.list$realisations <- realisations
	nrep <- temp.list$nrep <- length(levels(realisations))
	ind.rep <- split(1:n, realisations)
	vecdist <- function(x){as.vector(dist(x))}
	##
	## Initial values for parameters
	##
	## have to consider transformation, residuals from trend etc
#  var.data <- mean(tapply(data, realisations, var))
#  d.max <- max(by(ap$coords, ap$realisations, function(x) max(dist(x))))
#  if(missing(ini.cov.pars))
#    ini.cov.pars <- expand.grid(var.data/2, 3*var.data/4, var.data)
	if(any(class(ini.cov.pars) == "eyefit")){
		init <- nugget <- kappa <- NULL
		for(i in 1:length(ini.cov.pars)){
			init <- drop(rbind(init, ini.cov.pars[[i]]$cov.pars))
			nugget <- c(nugget, ini.cov.pars[[i]]$nugget)
			if(cov.model == "gneiting.matern")
				kappa <- drop(rbind(kappa, ini.cov.pars[[i]]$kappa))
			else
				kappa <- c(kappa, ini.cov.pars[[i]]$kappa)
		}
		ini.cov.pars <- init
	}
	if(any(class(ini.cov.pars) == "variomodel")){
		nugget <- ini.cov.pars$nugget
		kappa <- ini.cov.pars$kappa
		ini.cov.pars <- ini.cov.pars$cov.pars
	}
	if(is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)){
		ini.cov.pars <- as.matrix(ini.cov.pars)
		if(nrow(ini.cov.pars) == 1)
			ini.cov.pars <- as.vector(ini.cov.pars)
		else{
			if((cov.model != "pure.nugget") & (ncol(ini.cov.pars) != 2))
				stop("\nini.cov.pars must be a matrix or data.frame with 2 components: \ninitial values for sigmasq and phi")
		}
	}
	if(is.vector(ini.cov.pars)){
		if((cov.model != "pure.nugget") & (length(ini.cov.pars) != 2))
			stop("\nini.cov.pars must be a vector with 2 components: \ninitial values for sigmasq and phi")
	}
	##
	## Checking for multiple initial values for preliminar search of   
	## best initial value
	##
	thePos=-1
	
	if(is.matrix(ini.cov.pars) | (length(nugget) > 1) | (length(kappa) > 1) | (length(lambda) > 1) | (length(psiR) > 1) | (length(psiA) > 1)){
		if(messages.screen) cat("likfit: searching for best initial value ...")
		ini.temp <- matrix(ini.cov.pars, ncol=2)
		grid.ini <- as.matrix(expand.grid(sigmasq=unique(ini.temp[,1]), phi=unique(ini.temp[,2]), tausq=unique(nugget), kappa=unique(kappa), lambda=unique(lambda), psiR=unique(psiR), psiA=unique(psiA)))
		.likGRF.dists.vec= lapply(split(as.data.frame(coords), 
								realisations), vecdist)
				
		temp.f <- function(parms, coords, data, temp.list)
			return(loglik.GRF(geodata = geodata,
							coords = coords, data = as.vector(data),
							cov.model=temp.list$cov.model,
							cov.pars=parms[1:2],
							nugget=parms["tausq"], kappa=parms["kappa"],
							lambda=parms["lambda"], psiR=parms["psiR"],
							psiA=parms["psiA"], trend= trend,
							method.lik=temp.list$method.lik,
							compute.dists=FALSE,
							realisations = realisations))
		grid.lik <- apply(grid.ini, 1, temp.f, coords = coords,
				data = data, temp.list = temp.list)
		grid.ini <- grid.ini[(grid.lik != Inf) & (grid.lik != -Inf) & !is.na(grid.lik) & !is.nan(grid.lik),, drop=FALSE] 
		grid.lik <- grid.lik[(grid.lik != Inf) & (grid.lik != -Inf) & !is.na(grid.lik) & !is.nan(grid.lik)] 
		ini.temp <- grid.ini[which(grid.lik == max(grid.lik)),, drop=FALSE]
		if(all(ini.temp[,"phi"] == 0)) ini.temp <- ini.temp[1,, drop=FALSE]
		rownames(ini.temp) <- "initial.value"
		if(messages.screen){
			cat(" selected values:\n")
			print(rbind(format(ini.temp, digits=2), status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa, fix.lambda, fix.psiR, fix.psiA), "fix", "est")))
			cat(paste("likelihood value:", max(grid.lik), "\n"))
		}
		dimnames(ini.temp) <- NULL
		ini.cov.pars <- ini.temp[1:2]
		nugget <- ini.temp[3]
		kappa <- ini.temp[4]
		lambda <- ini.temp[5]
		psiR <- ini.temp[6]
		psiA <- ini.temp[7]
		grid.ini <- NULL
		.likGRF.dists.vec=NULL
	}
	##
	tausq <- nugget
	##
	## Box-Cox transformation for fixed lambda
	##
	if(fix.lambda) {
		if(abs(lambda - 1) < 0.0001) {
			temp.list$log.jacobian <- 0
			temp.list$z <- as.vector(data)
		} else {
			if(any(data <= 0))
				stop("Transformation option not allowed when there are zeros or negative data")
			Jdata <- data^(lambda - 1)
			if(any(Jdata <= 0))
				temp.list$log.jacobian <- log(prod(Jdata))
			else temp.list$log.jacobian <- sum(log(Jdata))
			Jdata <- NULL
			if(abs(lambda) < 0.0001)
				temp.list$z <- log(data)
			else temp.list$z <- ((data^lambda) - 1)/lambda
		}
	}else{
		temp.list$z <- as.vector(data)
		temp.list$log.jacobian <- NULL
	}
	##
	## Coordinates transformation for fixed anisotropy parameters
	##
	if(fix.psiR & fix.psiA){
		if(psiR != 1 | psiA != 0)
			coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
		 .likGRF.dists.vec= lapply(split(as.data.frame(coords), 
								realisations), vecdist) 
		range.dist <- range( .likGRF.dists.vec )
		max.dist <- max(range.dist)
		min.dist <- min(range.dist)
	}
	##
	##
	##
	ini <- ini.cov.pars[2]
	##  fixed.pars <- NULL
	lower.optim <- c(limits$phi["lower"])
	upper.optim <- c(limits$phi["upper"])
	fixed.values <- list()
	if(fix.nugget) {
		##    fixed.pars <- c(fixed.pars, 0)
		fixed.values$tausq <- nugget
	}	else {
		ini <- c(ini, nugget/ini.cov.pars[1])
		lower.optim <- c(lower.optim, limits$tausq.rel["lower"])
		upper.optim <- c(upper.optim, limits$tausq.rel["upper"])
	}
	if(fix.kappa){
		##    fixed.kappa <- c(fixed.pars, kappa)
		fixed.values$kappa <- kappa
	}	else {
		ini <- c(ini, kappa)
		lower.optim <- c(lower.optim, limits$kappa["lower"])
		upper.optim <- c(upper.optim, limits$kappa["upper"])
	}
	if(fix.lambda){
		##    fixed.pars <- c(fixed.pars, lambda)
		fixed.values$lambda <- lambda
	}	else {
		ini <- c(ini, lambda)
		lower.optim <- c(lower.optim, limits$lambda["lower"])
		upper.optim <- c(upper.optim, limits$lambda["upper"])
	}
	if(fix.psiR){
		##    fixed.pars <- c(fixed.pars, psiR)
		fixed.values$psiR <- psiR
	}	else {
		ini <- c(ini, psiR)
		lower.optim <- c(lower.optim, limits$psiR["lower"])
		upper.optim <- c(upper.optim, limits$psiR["upper"])
	}
	if(fix.psiA){
		##    fixed.pars <- c(fixed.pars, psiA)
		fixed.values$psiA <- psiA
	}	else {
		ini <- c(ini, psiA)
		lower.optim <- c(lower.optim, limits$psiA["lower"])
		upper.optim <- c(upper.optim, limits$psiA["upper"])
	}
	## This must be here, after the previous ones:
	if(fix.nugget & nugget > 0){
		## Warning: Inverting order here, ini will be now: c(phi,sigmasg)
		ini <- c(ini, ini.cov.pars[1])
		lower.optim <- c(lower.optim, limits$sigmasq["lower"])
		upper.optim <- c(upper.optim, limits$sigmasq["upper"])
		##    fixed.pars <- c(fixed.pars, ini.cov.pars[1])
		##    fixed.values$sigmasq <- 0
	}
	##
	names(ini) <- NULL
	if(length(ini) == 1) justone <- TRUE	else justone <- FALSE
	##
	ip <- list(f.tausq = fix.nugget, f.kappa = fix.kappa,
			f.lambda = fix.lambda,
			f.psiR = fix.psiR, f.psiA = fix.psiA)
	##
	npars <- beta.size + 2 + sum(unlist(ip)==FALSE)
	temp.list$coords <- coords
	temp.list$xmat <- split(as.data.frame(unclass(xmat)), realisations)
	temp.list$xmat <- lapply(temp.list$xmat, as.matrix)
	temp.list$n <- as.vector(unlist(lapply(temp.list$xmat, nrow)))
	##
	## Constant term in the likelihood
	##
	temp.list$loglik.cte <- rep(0, nrep)
	for(i in 1:nrep){
		if(method.lik == "ML"){
			if(ip$f.tausq & (tausq > 0))
				temp.list$loglik.cte[i] <-  (temp.list$n[i]/2)*(-log(2*pi))
			else
				temp.list$loglik.cte[i] <-  (temp.list$n[i]/2)*(-log(2*pi) +
							log(temp.list$n[i]) -1)
		}
		if(method.lik == "RML"){
			xx.eigen <- eigen(crossprod(temp.list$xmat[[i]]),
					symmetric = TRUE, only.values = TRUE)
			if(ip$f.tausq & (tausq > 0))
				temp.list$loglik.cte[i] <- - ((temp.list$n[i]-beta.size)/2)*(log(2*pi)) +
						0.5 * sum(log(xx.eigen$values))
			else
				temp.list$loglik.cte[i] <-  - ((temp.list$n[i]-beta.size)/2)*(log(2*pi)) +
						((temp.list$n[i]-beta.size)/2)*(log(temp.list$n[i]-beta.size)) -
						((temp.list$n[i]-beta.size)/2) + 0.5 * sum(log(xx.eigen$values))
		}
	}
		## Numerical minimization of the -loglikelihood
	##
	if(length(ini) == 1){
		if(upper.optim == Inf) upper.optim <- 50*max.dist
		# get rid of control arguments
		ldots = ldots[names(ldots)!= "control"]
		lik.minim <- do.call("optimize", c(list(.negloglik.GRF,
								lower=lower.optim,
								upper=upper.optim,
								fp=fixed.values,
								ip=ip, temp.list = temp.list,
								likGRF.dists.vec=
										.likGRF.dists.vec), 
						ldots))
		lik.minim <- list(par = lik.minim$minimum,
				value = lik.minim$objective,
				convergence = 0,
				message = "function optimize used")      
	}
	else{
		MET <- pmatch(names(ldots), names(formals(optim)))
		if(is.na(MET) || all(names(formals(optim))[MET] != "method"))
			ldots$method <- "L-BFGS-B"
		if(!is.null(names(ldots))){
			names(ldots)[which(as.logical(pmatch(names(ldots), "method", nomatch=0)))] <- "method"
		}
		if(!is.null(ldots$method) && ldots$method == "L-BFGS-B"){
			ldots$lower <- lower.optim
			ldots$upper <- upper.optim
		}
		
	
		if(!is.null(ldots$control$parscale)) {
			
			ldots$control$parscale = 
					c(ldots$control$parscale,
					rep(1, length(lower.optim)-
									length(ldots$control$parscale)))			
		}
		
 
	
		lik.minim <- do.call("optim", 
				c(list(par = ini, fn = .negloglik.GRF,
								fp=fixed.values, ip=ip, 
								temp.list = temp.list,
								likGRF.dists.vec=
										.likGRF.dists.vec), 
						ldots))
		##      lik.minim <- optim(par = ini, fn = .negloglik.GRF, method=optim.METHOD
		##                         lower=lower.optim, upper=upper.optim,
		##                         fp=fixed.values, ip=ip, temp.list = temp.list, ...)
	}
	lik.minim$control = ldots$control
	##
	if(messages.screen) cat("likfit: end of numerical maximisation.\n")
	par.est <- lik.minim$par
	if(any(par.est < 0)) par.est <- round(par.est, digits=12)
	phi <- par.est[1]
	##
	## Values of the maximised likelihood
	##
	if(is.R())
		loglik.max <- - lik.minim$value
	else
		loglik.max <- - lik.minim$objective
	##
	## Assigning values for estimated parameters
	##
	if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		psiA <- par.est[2]
	}
	if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		psiR <- par.est[2]
	}
	if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		psiR <- par.est[2]
		psiA <- par.est[3]
	}
	if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
		lambda  <- par.est[2]
	}
	if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		lambda  <- par.est[2]
		psiA <- par.est[3]
	}
	if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		lambda  <- par.est[2]
		psiR <- par.est[3]
	}
	if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		lambda  <- par.est[2]
		psiR <- par.est[3]
		psiA <- par.est[4]
	}
	if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
		kappa  <-  par.est[2]
	}
	if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		kappa  <-  par.est[2]
		psiA <- par.est[3]
	}
	if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		kappa  <-  par.est[2]
		psiR <- par.est[3]
	}
	if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		kappa  <-  par.est[2]
		psiR <- par.est[3]
		psiA <- par.est[4]
	}
	if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
		kappa <-  par.est[2]
		lambda <- par.est[3]
	}
	if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		kappa <-  par.est[2]
		lambda <- par.est[3]
		psiA <- par.est[4]
	}
	if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		kappa <-  par.est[2]
		lambda <- par.est[3]
		psiR<- par.est[4]
	}
	if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		kappa <-  par.est[2]
		lambda <- par.est[3]
		psiR<- par.est[4]
		psiA<- par.est[5]
	}
	if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- par.est[2]
	}
	if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- par.est[2]
		psiA<- par.est[3]
	}
	if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- par.est[2]
		psiR<- par.est[3]
	}
	if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- par.est[2]
		psiR<- par.est[3]
		psiA<- par.est[4]
	}
	if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- par.est[2]
		lambda <- par.est[3]
	}
	if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- par.est[2]
		lambda <- par.est[3]
		psiA <- par.est[4]
	}
	if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- par.est[2]
		lambda <- par.est[3]
		psiR <- par.est[4]
	}
	if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- par.est[2]
		lambda <- par.est[3]
		psiR <- par.est[4]
		psiA <- par.est[5]
	}
	if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- par.est[2]
		kappa <-  par.est[3]
	}
	if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- par.est[2]
		kappa <-  par.est[3]
		psiA <- par.est[4]
	}
	if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- par.est[2]
		kappa <-  par.est[3]
		psiR <- par.est[4]
	}
	if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- par.est[2]
		kappa <-  par.est[3]
		psiR <- par.est[4]
		psiA <- par.est[5]
	}
	if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- par.est[2]
		kappa <-  par.est[3]
		lambda <- par.est[4]
	}
	if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- par.est[2]
		kappa <-  par.est[3]
		lambda <- par.est[4]
		psiA <- par.est[5]
	}
	if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- par.est[2]
		kappa <-  par.est[3]
		lambda <- par.est[4]
		psiR <- par.est[5]
	}
	if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- par.est[2]
		kappa <-  par.est[3]
		lambda <- par.est[4]
		psiR <- par.est[5]
		psiA <- par.est[6]
	}
	##
	if(fix.nugget & nugget > 0){
		sigmasq <- par.est[length(par.est)]
		if(sigmasq > 1e-12) tausq <- nugget/sigmasq
		check.sigmasq <- TRUE
	}
	else check.sigmasq <- FALSE
	##
	##
	## Transforming data according to the estimated lambda (Box-Cox) parameter
	##
	if(!fix.lambda) {
		if(abs(lambda - 1) < 0.0001) {
			log.jacobian.max <- 0
		}
		else {
			if(any(data^(lambda - 1) <= 0))
				log.jacobian.max <- log(prod(data^(lambda - 1)))
			else log.jacobian.max <- sum(log(data^(lambda - 1)))
			temp.list$z <- ((data^lambda)-1)/lambda
		}
	}
	else{
		log.jacobian.max <- temp.list$log.jacobian
	}
	data.rep <- split(temp.list$z, realisations)
	coords.rep <- split(as.data.frame(coords), realisations)
	coords.rep <- lapply(coords.rep, as.matrix)
	##
	## Transforming coords for estimated anisotropy (if the case)
	##
	if(fix.psiR & fix.psiA)
		 .likGRF.dists.vec=NULL
	else{
		if(round(psiR, digits=6) != 1 | round(psiA, digits=6) != 0)
			coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
		rangevecdist <- function(x){range(as.vector(dist(x)))}
		range.dist <- lapply(split(as.data.frame(coords), realisations), rangevecdist)
		range.dist <- range(as.vector(unlist(range.dist)))
		max.dist <- max(range.dist)
		min.dist <- min(range.dist)
	}      
#  gc(verbose=FALSE)
	##
	## Computing estimated beta and tausq/sigmasq (if the case)
	##
	xivx <- matrix(0, ncol=beta.size, nrow=beta.size)
	xivy <- matrix(0, ncol=1, nrow=beta.size)
	yivy <- 0
	for(i in 1:nrep){
		ni <- temp.list$n[i]
		if((phi < 1e-12))
			V <- diag(x=(1+tausq), ni)
		else{
			if(check.sigmasq){
				if(sigmasq < 1e-12){
					if(!fix.nugget)
						V <- diag(x=(1+tausq), ni)
					else
						V <- diag(x=sqrt(tausq), ni)          
				}
				else
					V <- varcov.spatial(coords = coords.rep[[i]],
							cov.model = cov.model,
							kappa = kappa, nugget = tausq,
							cov.pars = c(1, phi))$varcov
			}
			else
				V <- varcov.spatial(coords = coords.rep[[i]],
						cov.model = cov.model,
						kappa = kappa, nugget = tausq,
						cov.pars = c(1, phi))$varcov
		}
		ivyx <- solve(V,cbind(data.rep[[i]],temp.list$xmat[[i]]))
		xivx <- xivx + crossprod(ivyx[,-1],temp.list$xmat[[i]])
		xivy <- xivy + crossprod(ivyx[,-1],data.rep[[i]])
		yivy <- yivy + crossprod(data.rep[[i]],ivyx[,1])
	}
	betahat <- geoR:::.solve.geoR(xivx, xivy)
	res <- as.vector(temp.list$z - xmat %*% betahat)
	if(!fix.nugget | (nugget < 1e-12)){
		ssres <- as.vector(yivy - 2*crossprod(betahat,xivy) +
						crossprod(betahat,xivx) %*% betahat)  
		if(method.lik == "ML")
			sigmasq <- ssres/n
		else
			sigmasq <- ssres/(n - beta.size)
	}
	if(fix.nugget){
		if(nugget > 0)
			tausq <- nugget
	}
	else tausq <- tausq * sigmasq
	betahat.var <- geoR:::.solve.geoR(xivx)
	if(sigmasq > 1e-12) betahat.var <- sigmasq * betahat.var
#  if(!fix.nugget & phi < 1e-16){
#    tausq <- sigmasq + tausq
#    sigmasq <- 0
#  }
	##
	## Preparing output
	##
	if((phi < 0.001*min.dist)){
		tausq <- tausq + sigmasq
		sigmasq <- 0
	}
	if((sigmasq < 1e-12)) phi <- 0
	##
	n.model.pars <- beta.size + 7
	par.su <- data.frame(status=rep(-9,n.model.pars))
	ind.par.su <- c(rep(0, beta.size), ip$f.tausq, 0, 0, ip$f.kappa,
			ip$f.psiR, ip$f.psiA,ip$f.lambda)
	par.su$status <- ifelse(ind.par.su,"fixed", "estimated")
	par.su$values <- round(c(betahat, tausq, sigmasq, phi, kappa, psiR, psiA, lambda), digits=4)
	if(beta.size == 1) beta.name <- "beta"
	else beta.name <- paste("beta", 0:(beta.size-1), sep="")
	row.names(par.su) <- c(beta.name, "tausq", "sigmasq", "phi", "kappa",
			"psiR", "psiA", "lambda")
	par.su <- par.su[c((1:(n.model.pars-3)), n.model.pars-1, n.model.pars-2, n.model.pars),] 
	##
	lik.results <- list(cov.model = cov.model,
			nugget = tausq,
			cov.pars=c(sigmasq, phi),
			sigmasq = sigmasq,
			phi = phi,
			kappa = kappa,
			beta = as.vector(betahat),
			beta.var = betahat.var,
			lambda = lambda,
			aniso.pars = c(psiA = psiA, psiR = psiR),
			tausq = tausq,
			practicalRange = practicalRange(cov.model=cov.model,
					phi = phi, kappa = kappa),
			method.lik = method.lik, trend = trend,
			loglik = loglik.max,
			npars = npars,
			AIC = -2 * (loglik.max - npars),
			BIC = -2 * (loglik.max - 0.5 * log(n) * npars),
#                      residuals = res,
			parameters.summary = par.su,
			info.minimisation.function = lik.minim,
			max.dist = max.dist,
			trend = trend,
			trend.matrix= xmat,
			transform.info = list(fix.lambda = fix.lambda,
					log.jacobian = log.jacobian.max))
	##
	## Likelihood results for the model without spatial correlation
	##
	if(nospatial){
		if(fix.lambda){
			beta.ns <- geoR:::.solve.geoR(crossprod(xmat), crossprod(xmat, temp.list$z))
			ss.ns <- sum((as.vector(temp.list$z - xmat %*% beta.ns))^2)
			if(method.lik == "ML"){
				nugget.ns <- ss.ns/n
				loglik.ns <- (n/2)*((-log(2*pi)) - log(nugget.ns) - 1) + temp.list$log.jacobian
			}
			if(method.lik == "RML"){
				nugget.ns <- ss.ns/(n-beta.size)
				loglik.ns <- ((n-beta.size)/2)*((-log(2*pi)) - log(nugget.ns) -1) +
						temp.list$log.jacobian
			}
			npars.ns <- beta.size + 1 + !fix.lambda
			lambda.ns <- lambda
		}
		else{
			warning("this part hasn't been written in geostatsp, call geoR directly instead")
			if(is.R())
				lik.lambda.ns <- optim(par=1, fn = geoR:::.negloglik.boxcox,
						method = "L-BFGS-B",
						lower = limits$lambda["lower"],
						upper = limits$lambda["upper"],
						data = data, xmat = xmat,
						lik.method = method.lik)
			else
				lik.lambda.ns <- nlminb(par=1, fn = geoR:::.negloglik.boxcox,
						lower=limits$lambda["lower"],
						upper=limits$lambda["upper"],
						data = data, xmat = xmat,
						lik.method = method.lik)
			lambda.ns <- lik.lambda.ns$par
			if(abs(lambda) < 0.0001) tdata.ns <- log(data)
			else tdata.ns <- ((data^lambda.ns)-1)/lambda.ns
			beta.ns <- geoR:::.solve.geoR(crossprod(xmat),crossprod(xmat,tdata.ns))
			ss.ns <- sum((as.vector(tdata.ns - xmat %*% beta.ns))^2)
			if(is.R())
				value.min.ns <- lik.lambda.ns$value
			else
				value.min.ns <- lik.lambda.ns$objective
			if(method.lik == "ML"){
				loglik.ns <- (- value.min.ns)+ (n/2)*((-log(2*pi)) + log(n) - 1)
				nugget.ns <- ss.ns/n
			}
			if(method.lik == "RML"){
				nugget.ns <- ss.ns/(n-beta.size)
				loglik.ns <- (- value.min.ns)+ ((n-beta.size)/2)*((-log(2*pi)) +
							log(n-beta.size) - 1)
			}      
			npars.ns <- beta.size + 1 + !fix.lambda
		}
		lik.results$nospatial <- list(beta.ns = beta.ns, variance.ns = nugget.ns,
				loglik.ns = loglik.ns, npars.ns = npars.ns,
				lambda.ns = lambda.ns, AIC.ns = -2 * (loglik.ns - npars.ns),
				BIC.ns = -2 * (loglik.ns - 0.5 * log(n) * npars.ns))
	}
	##
	## Assigning names to the components of the mean vector beta
	##
	if(length(lik.results$beta.var) == 1)
		lik.results$beta.var <- as.vector(lik.results$beta.var)
	if(length(lik.results$beta) > 1){
		##    if(inherits(trend, "formula") || (!is.null(class(trend)) && any(class(trend) == "trend.spatial")))
		if(inherits(trend, "formula") || (length(class(trend)) > 0 && any(class(trend) == "trend.spatial")))
			beta.names <- c("intercept", paste("covar", 1:(ncol(xmat)-1), sep = ""))
		else
		if(trend == "1st")
			beta.names <- c("intercept", "x", "y")
		else
		if(trend == "2nd")
			beta.names <- c("intercept", "x", "y", "x2", "xy", "y2")
		names(lik.results$beta) <- beta.names
	}
	##
	## Computing residuals and predicted values
	## (isolated components of the model)
	##
	if(components) {
		if(!fix.psiR & !fix.psiA)
			if(psiR != 1 | psiA != 0)
				coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
		#coords.rep <- split(as.data.frame(coords), realisations)
		#res.rep <- split(res, realisations)
		trend.comp <- temp.list$z - res
		spatial.comp <- list()
		for(i in 1:nrep){
#      invcov <- varcov.spatial(coords = coords[ind.rep[[i]],], cov.model = cov.model, 
#                               kappa = kappa, nugget = tausq,
#                               cov.pars = c(sigmasq, phi), inv=TRUE)$inverse 
#      covmat.signal <- varcov.spatial(coords = coords[ind.rep[[i]],],
#                                      cov.model = cov.model, 
#                                      kappa = kappa, nugget = 0,
#                                      cov.pars = c(sigmasq, phi))$varcov
			spatial.comp[[i]] <- as.vector(varcov.spatial(coords = coords[ind.rep[[i]],],
									cov.model = cov.model, 
									kappa = kappa, nugget = 0,
									cov.pars = c(sigmasq, phi))$varcov %*%
							varcov.spatial(coords = coords[ind.rep[[i]],],
									cov.model = cov.model, 
									kappa = kappa, nugget = tausq,
									cov.pars = c(sigmasq, phi), inv=TRUE)$inverse %*%
							res[ind.rep[[i]]]) 
		}
		spatial.comp <- as.vector(unlist(spatial.comp))[as.vector(unlist(ind.rep))]
		predict.comp <- trend.comp + spatial.comp
		residual.comp <- as.vector(temp.list$z - predict.comp)
#    residual.std <- as.vector(invcov %*% residual.comp)
#    residual.trend.std <- as.vector(invcov %*% res)
		lik.results$model.components <-
				data.frame(trend = trend.comp, spatial = spatial.comp, residuals = residual.comp)
#    lik.results$s2.random <- (crossprod(res,invcov) %*% res)/(n - beta.size)
#    lik.results$s2 <- (crossprod(residual.comp,invcov) %*% residual.comp)/(n - beta.size)
	}
	##
	lik.results$contrasts <- xmat.contrasts
	lik.results$call <- call.fc
	##
	## Assigning classes
	##
	oldClass(lik.results) <- c("likGRF", "variomodel")
	##
	## Some warning messages about particular possible results
	##
	if(messages.screen){
		if((lik.results$cov.pars[1] < (0.01 * (lik.results$nugget + lik.results$cov.pars[1])))& lik.results$cov.pars[2] > 0)
			cat("\nWARNING: estimated sill is less than 1 hundredth of the total variance. Consider re-examine the model excluding spatial dependence\n" )      
		if((lik.results$cov.pars[2] > (10 * max.dist)) & lik.results$cov.pars[1] > 0 )
			cat("\nWARNING: estimated range is more than 10 times bigger than the biggest distance between two points. Consider re-examine the model:\n 1) excluding spatial dependence if estimated sill is too low and/or \n 2) taking trends (covariates) into account\n" ) 
		if(((lik.results$cov.pars[2] < (0.1 * min.dist)) & (lik.results$cov.pars[1] > 0)) & lik.results$cov.pars[2] > 0)
			cat("\nWARNING: estimated range is less than 1 tenth of the minimum distance between two points. Consider re-examine the model excluding spatial dependence\n" ) 
	}
	##
	attr(lik.results, "geodata") <- name.geodata

	# add names of covariates
	names(lik.results$beta) = colnames(xmat)
	if(is.matrix(lik.results$beta.var)){
		dimnames(lik.results$beta.var) = list(colnames(xmat),colnames(xmat))

		beta.table = cbind(
			est = lik.results$beta,
			std.err = sqrt(diag(lik.results$beta.var))
			)
	 } else {
		 beta.table = cbind(
				 est = lik.results$beta,
				 std.err = sqrt(lik.results$beta.var)
		 )
	 }
		
	beta.table= cbind(beta.table,
			z = beta.table[,"est"]/beta.table[,"std.err"])
	beta.table= cbind(beta.table,
		p = pnorm(abs(beta.table[,"z"]), lower.tail=FALSE)		
	)
	beta.table = rbind(beta.table, range = c(
					lik.results$phi * sqrt(8*lik.results$kappa),
					rep(NA, dim(beta.table)[2]-1)
					))
	lik.results$beta.table=beta.table
	return(lik.results)
}

