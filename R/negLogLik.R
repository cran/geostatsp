
".negloglik.GRF" <-
		function(pars, fp, ip, temp.list, likGRF.dists.vec)
	# this function is a modified version of the function of the same name 
	# in the geoR package by Paulo Ribeiro and Peter Diggle
	
	
	### pars : values for the parameters to be estimated
## sequence is c(phi, tausq, kappa, lambda, psiR, psiA, sigmasq)
### fixed pars: parameters considered fixed
### ind.pars : list indicating which are fixed and which are to be estimated
##
## Warning:
##  if fix.nugget = TRUE and nugget > 0 ,
## sigmasq should be passed and fp$nugget is the value of the nugget
## otherwise the RELATIVE nugget should be passed
{
	p <- temp.list$beta.size
	log.jacobian <- temp.list$log.jacobian
	## Obligatory parameter:
	phi <- pars[1]
	## Others
	if(ip$f.tausq){
		if(fp$tausq > 0){
			npars.min <- length(pars)
			sigmasq <- pars[npars.min]
		}
		else sigmasq <- 1
	}
	else sigmasq <- 1
	if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- fp$tausq
		kappa <- fp$kappa
		lambda <- fp$lambda
		psiR <- fp$psiR
		psiA <- fp$psiA
	}
	if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- fp$tausq
		kappa <- fp$kappa
		lambda <- fp$lambda
		psiR <- fp$psiR
		psiA <- pars[2]
	}
	if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- fp$tausq
		kappa <- fp$kappa
		lambda <- fp$lambda
		psiR <- pars[2]
		psiA <- fp$psiA
	}
	if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- fp$tausq
		kappa <- fp$kappa
		lambda <- fp$lambda
		psiR <- pars[2]
		psiA <- pars[3]
	}
	if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- fp$tausq
		kappa <- fp$kappa
		lambda <- pars[2]
		psiR <- fp$psiR
		psiA <- fp$psiA
	}
	if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- fp$tausq
		kappa <- fp$kappa
		lambda <- pars[2]
		psiR <- fp$psiR
		psiA <- pars[3]
	}
	if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- fp$tausq
		kappa <- fp$kappa
		lambda <- pars[2]
		psiR <- pars[3]
		psiA <- fp$psiA
	}
	if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- fp$tausq
		kappa <- fp$kappa
		lambda <- pars[2]
		psiR <- pars[3]
		psiA <- pars[4]
	}
	if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- fp$tausq
		kappa <- pars[2]
		lambda <- fp$lambda
		psiR <- fp$psiR
		psiA <- fp$psiA
	}
	if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- fp$tausq
		kappa <- pars[2]
		lambda <- fp$lambda
		psiR <- fp$psiR
		psiA <- pars[3]
	}
	if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- fp$tausq
		kappa <- pars[2]
		lambda <- fp$lambda
		psiR <- pars[3]
		psiA <- fp$psiA
	}
	if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- fp$tausq
		kappa <- pars[2]
		lambda <- fp$lambda
		psiR <- pars[3]
		psiA <- pars[4]
	}
	if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- fp$tausq
		kappa <- pars[2]
		lambda <- pars[3]
		psiR <- fp$psiR
		psiA <- fp$psiA
	}
	if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- fp$tausq
		kappa <- pars[2]
		lambda <- pars[3]
		psiR <- fp$psiR
		psiA <- pars[4]
	}
	if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- fp$tausq
		kappa <- pars[2]
		lambda <- pars[3]
		psiR <- pars[4]
		psiA <- fp$psiA
	}
	if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- fp$tausq
		kappa <- pars[2]
		lambda <- pars[3]
		psiR <- pars[4]
		psiA <- pars[5]
	}
	if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- pars[2]
		kappa <- fp$kappa
		lambda <- fp$lambda
		psiR <- fp$psiR
		psiA <- fp$psiA
	}
	if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- pars[2]
		kappa <- fp$kappa
		lambda <- fp$lambda
		psiR <- fp$psiR
		psiA <- pars[3]
	}
	if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- pars[2]
		kappa <- fp$kappa
		lambda <- fp$lambda
		psiR <- pars[3]
		psiA <- fp$psiA
	}
	if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- pars[2]
		kappa <- fp$kappa
		lambda <- fp$lambda
		psiR <- pars[3]
		psiA <- pars[4]
	}
	if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- pars[2]
		kappa <- fp$kappa
		lambda <- pars[3]
		psiR <- fp$psiR
		psiA <- fp$psiA
	}
	if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- pars[2]
		kappa <- fp$kappa
		lambda <- pars[3]
		psiR <- fp$psiR
		psiA <- pars[4]
	}
	if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- pars[2]
		kappa <- fp$kappa
		lambda <- pars[3]
		psiR <- pars[4]
		psiA <- fp$psiA
	}
	if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- pars[2]
		kappa <- fp$kappa
		lambda <- pars[3]
		psiR <- pars[4]
		psiA <- pars[5]
	}
	if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- pars[2]
		kappa <- pars[3]
		lambda <- fp$lambda
		psiR <- fp$psiR
		psiA <- fp$psiA
	}
	if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- pars[2]
		kappa <- pars[3]
		lambda <- fp$lambda
		psiR <- fp$psiR
		psiA <- pars[4]
	}
	if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- pars[2]
		kappa <- pars[3]
		lambda <- fp$lambda
		psiR <- pars[4]
		psiA <- fp$psiA
	}
	if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- pars[2]
		kappa <- pars[3]
		lambda <- fp$lambda
		psiR <- pars[4]
		psiA <- pars[5]
	}
	if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
		tausq <- pars[2]
		kappa <- pars[3]
		lambda <- pars[4]
		psiR <- fp$psiR
		psiA <- fp$psiA
	}
	if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
		tausq <- pars[2]
		kappa <- pars[3]
		lambda <- pars[4]
		psiR <- fp$psiR
		psiA <- pars[5]
	}
	if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
		tausq <- pars[2]
		kappa <- pars[3]
		lambda <- pars[4]
		psiR <- pars[5]
		psiA <- fp$psiA
	}
	if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
		tausq <- pars[2]
		kappa <- pars[3]
		lambda <- pars[4]
		psiR <- pars[5]
		psiA <- pars[6]
	}
	##
	if(temp.list$print.pars){
		running.pars <- c(phi = phi, tausq = tausq, kappa =kappa, psiA = psiA, psiR = psiR, lambda = lambda)
		if(ip$f.tausq && fp$tausq > 0)
			running.pars <- c(sigmasq=sigmasq, running.pars)
		print(running.pars)
	}
	##
	## Absurd values
	##
	if(kappa < 1e-04 | (tausq+sigmasq) < (.Machine$double.eps^0.5) |
			any(c(phi, tausq, sigmasq, kappa) < 0))
		return(.Machine$double.xmax^0.5)
	##
	## Anisotropy
	##
	if(!ip$f.psiR | !ip$f.psiA){
		coords.c <- coords.aniso(temp.list$coords, aniso.pars=c(psiA, psiR))
		vecdist <- function(x){as.vector(dist(x))}
		assign("likGRF.dists.vec", lapply(split(as.data.frame(coords.c),
								temp.list$realisations), vecdist))
	}
	##
	## Box-Cox transformation
	##
	if(!ip$f.lambda){
		if(abs(lambda - 1) < 0.0001) {
			log.jacobian <- 0
		}
		else {
			if(any(temp.list$z <= 0))
				stop("Transformation not allowed for zero or negative data")
			data <- temp.list$z^(lambda - 1)
			if(any(data <= 0)) log.jacobian <- log(prod(data))
			else log.jacobian <- sum(log(data))
			data <- NULL
		}
		if(abs(lambda) < 0.0001)
			data <- log(temp.list$z)
		else data <- ((temp.list$z^lambda) - 1)/lambda
	}
	else data <- temp.list$z
	data <- split(data, as.factor(temp.list$realisations))
	##
	## Computing likelihood
	##
	sumnegloglik <- 0
	for(i in 1:temp.list$nrep){
		## NOTE: Likelihood for Independent observations 
		##       arbitrary criteria used here:
		##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
		##
		n <- temp.list$n[i]
		xmat <- temp.list$xmat[[i]]
		z <- data[[i]]
		if((phi < 1e-16) | (sigmasq < 1e-16)){
			if(ip$f.tausq)
				v <- list(varcov = diag(x=(tausq+sigmasq), n),
						log.det.to.half = (n/2) * log(tausq+sigmasq))
			else
				v <- list(varcov = diag(x=(1+tausq), n),
						log.det.to.half = (n/2) * log(1+tausq))
		}
		else
			v <- varcov.spatial(dists.lowertri = get("likGRF.dists.vec")[[i]],
					cov.model = temp.list$cov.model, kappa=kappa,
					nugget = tausq, cov.pars=c(sigmasq, phi),
					det = TRUE)
		if(!is.null(v$crash.parms)) return(.Machine$double.xmax^0.5)
		ivx <- solve(v$varcov,xmat)
		xivx <- crossprod(ivx,xmat)
		betahat <- try(.solve.geoR(xivx,crossprod(ivx,z)), silent=TRUE)
		if(inherits(betahat, "try-error")){
			t.ei <- eigen(xivx, symmetric = TRUE)
#      if(exists("trySilent"))
			betahat <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% crossprod(ivx,z), silent=TRUE)
#      else{
#        error.now <- options()$show.error.message
#        options(show.error.messages = FALSE)
#        betahat <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% crossprod(ivx,z))
#        if(is.null(error.now) || error.now) options(show.error.messages = TRUE)        
#      }
		}
		if(inherits(betahat, "try-error"))
			stop("Covariates have very different orders of magnitude. Try to multiply and/or divide them to bring them to similar orders of magnitude") 
		res <- z - xmat %*% betahat
		ssres <- drop(crossprod(res,solve(v$varcov,res)))
		if(temp.list$method.lik == "ML"){
			if(ip$f.tausq & (tausq > 0))
				negloglik <- v$log.det.to.half +  0.5 * ssres
			else
				negloglik <- (n/2) * log(ssres) +  v$log.det.to.half
		}
		if(temp.list$method.lik == "RML"){
			if(length(as.vector(xivx)) == 1) {
				choldet <- 0.5 * log(xivx)
			}
			else {
				chol.xivx <- chol(xivx)
				choldet <- sum(log(diag(chol.xivx)))
			}
			if(ip$f.tausq & (tausq > 0))
				negloglik <- v$log.det.to.half +  0.5 * ssres + choldet
			else
				negloglik <- ((n-p)/2) * log(ssres) +  v$log.det.to.half + choldet
		}  
		negloglik <- negloglik - temp.list$loglik.cte[i]
		sumnegloglik <- sumnegloglik + negloglik
	}
	sumnegloglik <- sumnegloglik - log.jacobian
	if(sumnegloglik > (.Machine$double.xmax^0.5) | sumnegloglik == Inf | sumnegloglik == -Inf)
		sumnegloglik <- .Machine$double.xmax^0.5
	if(temp.list$print.pars)
		cat(paste("log-likelihood = ", -sumnegloglik, "\n"))
	return(sumnegloglik) 
}
".solve.geoR" <-
		function (a, b = NULL, ...) 
{
	# this function is a modified version of the function of the same name 
	# in the geoR package by Paulo Ribeiro and Peter Diggle
	
	
	a <- eval(a)
	b <- eval(b)
#  if(exists("trySilent")){
	if (is.null(b)) res <- try(solve(a, ...), silent=TRUE)
	else res <- try(solve(a, b, ...), silent=TRUE)
#  }
#  else{
#    error.now <- options()$show.error.messages
#    if (is.null(error.now) | error.now) 
#      on.exit(options(show.error.messages = TRUE))
#    options(show.error.messages = FALSE)
#    if (is.null(b)) res <- try(solve(a, ...))
#    else res <- try(solve(a, b, ...))
#  }
	if (inherits(res, "try-error")) {
		test <- all.equal.numeric(a, t(a), 100 * .Machine$double.eps)
		if(!(is.logical(test) && test)){
			##      options(show.error.messages = TRUE)
			stop("matrix `a' is not symmetric")
		}
		t.ei <- eigen(a, symmetric = TRUE)
#    if(exists("trySilent")){
		if (is.null(b))
			res <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)), silent=TRUE)
		else
			res <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% b, silent=TRUE)
#    }
#    else{
#      if (is.null(b)) res <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)))
#      else res <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% b)
#    }
		if (any(is.na(res)) | any(is.nan(res)) | any(is.infinite(res))) 
			oldClass(res) <- "try-error"
	}
	if (inherits(res, "try-error")) 
		stop("Singular matrix. Covariates may have different orders of magnitude.")
	return(res)
}

".geoR.cov.models" <-
		c("matern", "exponential", "gaussian", "spherical",
				"circular", "cubic", "wave", "linear", "power",
				"powered.exponential", "stable", "cauchy", "gencauchy",
				"gneiting", "gneiting.matern", "pure.nugget")

"geoRCovModels" <- .geoR.cov.models

".check.geoRparameters.values" <- function(list, messages = TRUE)
{
	if(!is.null(list$nugget))
		if(list$nugget < 0) stop("value for nugget must be non-negative")
	if(!is.null(list$psiA))
		if(list$psiA < 0) stop("value for psiA must be non-negative")
	if(!is.null(list$psiR))
		if(list$psiR < 1) stop("value for psiA must be >= 1")
	if(!is.null(list$kappa)){
		if(is.null(list$cov.model)) stop("cov.model must be provided when checking values of kappa")
		if(list$kappa[1] <= 0) stop("parameter kappa[1] must be greater than 0")
		if(any(list$cov.model == c("powered.exponential", "matern", "gneiting.matern", "gencauchy","cauchy"))){
			if(any(list$cov.model == c("gneiting.matern", "gencauchy"))){
				if(length(list$kappa) != 2)
					stop(paste("kappa must be of length 2 for the", list$cov.model, "correlation function"))
				if(list$cov.model == "gencauchy" && (list$kappa[2] <=0 | list$kappa[2] >2))
					stop("for the gencauchy model the kappa must be within (0,2]")          
			}
			else{
				if(list$cov.model == "powered.exponential" && (list$kappa <=0 | list$kappa >2))
					stop("for the powered.exponential model the kappa must be within (0,2]")          
				if(list$cov.model == "stable" && (list$kappa <=0 | list$kappa >2))
					stop("for the stable model the kappa must be within (0,2]")          
			}
		}
		else{
			if(messages && !is.null(list$kappa))
				cat(paste("kappa not used for the",list$cov.model, "correlation function\n"))
		}
	}
	return(invisible())
}
".negloglik.boxcox" <-
		function(lambda.val, data, xmat, lik.method = "ML")
{
	if(length(lambda.val) == 2){
		data <- data + lambda.val[2]
		lambda <- lambda.val[1]
	}
	else lambda <- lambda.val
	lambda <- unname(lambda)
	n <- length(data)
	beta.size <- ncol(xmat)
	if(isTRUE(all.equal(unname(lambda), 0))) yt <- log(data)
	else yt <- ((data^lambda) - 1)/lambda
	beta <- solve(crossprod(xmat), crossprod(xmat, yt))
	ss <- sum((drop(yt) - drop(xmat %*% beta))^2)
	if(lik.method == "ML")
		neglik <- (n/2) * log(ss) - ((lambda - 1) * sum(log(data)))
	if(lik.method == "RML"){
		xx <- crossprod(xmat)
		if(length(as.vector(xx)) == 1)
			choldet <- 0.5 * log(xx)
		else
			choldet <- sum(log(diag(chol(xx))))
		neglik <- ((n-beta.size)/2) * log(ss) + choldet -
				((lambda - 1) * sum(log(data)))
	}
	if(mode(neglik) != "numeric") neglik <- Inf
	return(drop(neglik))
}

