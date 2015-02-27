\name{likfitLgm}
\alias{likfitLgm}
\alias{loglikLgm}

\title{Likelihood Based Parameter Estimation for Gaussian Random Fields}

\description{
  \emph{Maximum likelihood} (ML) or \emph{restricted maximum likelihood} (REML)
  parameter estimation for (transformed) Gaussian random fields.
}



\usage{


likfitLgm(formula, data, 
		paramToEstimate = c("range","nugget"),
		reml=TRUE,
		coordinates=data,
		param=NULL,
		upper=NULL,lower=NULL, parscale=NULL,
		verbose=FALSE)

		
loglikLgm(param, 
		data, formula, coordinates=data,
		reml=TRUE, 
		minustwotimes=TRUE,
		moreParams=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{formula}{
A formula for the fixed effects portion of the model, specifying a response and covariates. 
Alternately, \code{data} can be a vector of observations and \code{formula} can be 
a model matrix.
}
 \item{data}{
	An object of class \code{SpatialPointsDataFrame}, a vector of observations, 
	or a data frame containing observations and covariates.
}
\item{coordinates}{A \code{SpatialPoints} object containing the locations of each observation, 
which defaults to \code{data}.  Alternately, \code{coordinates} can be a
 \code{\link[Matrix]{symmetricMatrix-class}} or \code{dist} object 
reflecting the distance matrix of these coordinates (though this is only permitted if the 
model is isotropic).}
\item{param}{A vector of model parameters, with named elements being amongst 
\code{range, nugget, boxcox, shape, anisoAngleDegrees, anisoAngleRadians,
anisoRatio}, and possibly \code{variance} (see \code{\link{matern}}). 
 When calling \code{likfitLgm} this vector is a 
combination of starting values for parameters to be estiamated and fixed values of parameters
which will not be estimated.  For \code{loglikLgm}, it is the covariance parameters 
for which the likelihood will be evaluated.
}
\item{reml}{Whether to use Restricted Likelihood rather than Likelihood, defaults to \code{TRUE}.}
\item{paramToEstimate}{Vector of names of model parameters to estimate, with parameters
excluded from this list being fixed.  The variance parameter and regression coefficients
are always estimated even if not listed.}
\item{lower}{Named vector of lower bounds for model parameters passed to \code{\link[stats]{optim}}, defaults are 
used for parameters not specified.}
\item{upper}{Upper bounds, as above.}
\item{parscale}{Named vector of scaling of parameters passed as \code{control=list(parscale=parscale)} 
to  \code{\link[stats]{optim}}.}
\item{minustwotimes}{Return -2 times the log likelihood rather than the likelihood}
\item{moreParams}{Vector of additional parameters, combined with \code{param}.  
Used for passing fixed parameters to \code{loglikLgm} from within \code{\link[stats]{optim}}.
}
\item{verbose}{if \code{TRUE} information is printed by \code{\link[stats]{optim}}.}
}

\value{
\code{likfitLgm} produces list with elements
\item{parameters}{Maximum Likelihood Estimates of model parameters}
\item{varBetaHat}{Variance matrix of the estimated regression parameters}
\item{optim}{results from \code{\link[stats]{optim}} }
\item{trend}{Either formula for the fixed effects or names of the columns
of the model matrix, depending on \code{trend} supplied.}
\item{summary}{a table of parameter estimates, standard errors,  confidence intervals, p values, and 
  	a logical value indicating whether each parameter was estimated as opposed to fixed.}
\item{resid}{residuals, being the observations minus the fixed effects, on the
transformed scale.}
\code{loglikLgm} returns a scalar value, either the log likelihood or -2 times the
log likelihood.  Attributes of this result include the vector of
parameters (including the MLE's computed for the variance and coefficients),
and the variance matrix of the coefficient MLE's.
 }


\seealso{
\code{\link{lgm}}
}
\examples{
n=40
mydat = SpatialPointsDataFrame(
	cbind(runif(n), seq(0,1,len=n)), 
	data=data.frame(cov1 = rnorm(n), cov2 = rpois(n, 0.5))
	)

# simulate a random field
trueParam = c(variance=2^2, range=0.35, shape=2, nugget=0.5^2)
set.seed(1)

mydat@data= cbind(mydat@data , RFsimulate(model=trueParam,x=mydat)@data)

# add fixed effects
mydat$Y = -3 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
	mydat$sim + rnorm(length(mydat), 0, sd=sqrt(trueParam["nugget"]))

spplot(mydat, "sim", col.regions=rainbow(10), main="U")
spplot(mydat, "Y", col.regions=rainbow(10), main="Y")


myres = likfitLgm(
	formula=Y ~ cov1 + cov2, 
	data=mydat, 
	param=c(range=0.1,nugget=0.1,shape=2), 
	paramToEstimate = c("range","nugget")
	)



myres$summary[,1:4]

if(interactive()  | Sys.info()['user'] =='patrick') {

# plot variograms of data, true model, and estimated model
myv = variog(mydat, formula=Y ~ cov1 + cov2,option="bin", max.dist=0.5)
# myv will be NULL if geoR isn't installed
if(!is.null(myv)){
plot(myv, ylim=c(0, max(c(1.2*sum(trueParam[c("variance", "nugget")]),myv$v))),
	main="variograms")
distseq = seq(0, 0.5, len=50)
lines(distseq, 
	sum(myres$param[c("variance", "nugget")]) - matern(distseq, param=myres$param),
	col='blue', lwd=3)
lines(distseq, 
	sum(trueParam[c("variance", "nugget")]) - matern(distseq, param=trueParam),
	col='red')	

legend("bottomright", fill=c("black","red","blue"), 
	legend=c("data","true","MLE"))
}

# without a nugget
myresNoN = likfitLgm(
	formula=Y ~ cov1 + cov2, 
	data=mydat, 
	param=c(range=0.1,nugget=0,shape=1), 
	paramToEstimate = c("range")
	)

myresNoN$summary[,1:4]


# plot variograms of data, true model, and estimated model
myv = variog(mydat, formula=Y ~ cov1 + cov2,option="bin", max.dist=0.5)

if(!is.null(myv)){
plot(myv, ylim=c(0, max(c(1.2*sum(trueParam[c("variance", "nugget")]),myv$v))),
	main="variograms")
	
distseq = seq(0, 0.5, len=50)
lines(distseq, 
	sum(myres$param[c("variance", "nugget")]) - matern(distseq, param=myres$param),
	col='blue', lwd=3)
lines(distseq, 
	sum(trueParam[c("variance", "nugget")]) - matern(distseq, param=trueParam),
	col='red')	

lines(distseq, 
	sum(myresNoN$param[c("variance", "nugget")]) - 
			matern(distseq, param=myresNoN$param),
	col='green', lty=2, lwd=3)	
legend("bottomright", fill=c("black","red","blue","green"), 
	legend=c("data","true","MLE","no N"))
}


# calculate likelihood
temp=loglikLgm(param=myres$param, 
		data=mydat, 
		formula = Y ~ cov1 + cov2,
		reml=FALSE, minustwotimes=FALSE)



# an anisotropic example


trueParamAniso = param=c(variance=2^2, range=0.2, shape=2,
	 nugget=0,anisoRatio=4,anisoAngleDegrees=10, nugget=0)

mydat$U = geostatsp::RFsimulate(trueParamAniso,mydat)$sim

mydat$Y = -3 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
	mydat$U + rnorm(length(mydat), 0, sd=sqrt(trueParamAniso["nugget"]))



par(mfrow=c(1,2))
oldmar = par("mar")
par(mar=rep(0.1, 4))

plot(mydat, col=as.character(cut(mydat$U, breaks=50, labels=heat.colors(50))),
	pch=16, main="aniso")
 
plot(mydat, col=as.character(cut(mydat$Y, breaks=50, labels=heat.colors(50))),
	pch=16,main="iso")
par(mar=oldmar)


myres = likfitLgm( 
	formula=Y ~ cov1 + cov2, 
	data=mydat,
	param=c(range=0.1,nugget=0,shape=2, anisoAngleDegrees=0, anisoRatio=2), 
	paramToEstimate = c("range","nugget","anisoRatio","anisoAngleDegrees") 
	)

myres$summary

par(mfrow=c(1,2))

myraster = raster(nrows=30,ncols=30,xmn=0,xmx=1,ymn=0,ymx=1)
covEst = matern(myraster, y=c(0.5, 0.5), par=myres$param)
covTrue = matern(myraster, y=c(0.5, 0.5), par=trueParamAniso)

plot(covEst, main="estimate")
plot(covTrue, main="true")

par(mfrow=c(1,1))
}
}

