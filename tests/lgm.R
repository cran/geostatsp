library("geostatsp")
data("swissRain")
library("rgdal")



bob = function(x) {
	thepar = x$param
	pdf(tempfile("lgm", tmpdir=".", fileext=".pdf"))
	plot(x$predict[["predict"]], main=
					paste(
							paste(names(thepar), thepar, sep="="),
							collapse=", "),cex.main=0.3
	)
	dev.off()
}



# specify formula name of raster layer
swissFit = lgm(data=swissRain, formula=rain~ CHE_alt,
		locations=80, covariates=swissAltitude,
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFit)
swissFit$param
bob(swissFit)


# specify formula using name of list element

swissFitAgain = lgm(data=swissRain, formula=rain~ elev+land,
		locations=80, covariates=list(elev=swissAltitude,land=swissLandType),
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param
bob(swissFitAgain)

swissFitAgain = lgm(data=swissRain, formula="rain",
		locations=80, covariates=swissAltitude,
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param
bob(swissFitAgain)


swissFitAgain = lgm(data=swissRain, formula="rain",
		locations=80, covariates=list(elev=swissAltitude,land=swissLandType),
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param
bob(swissFitAgain)


# land type, factor covariate
swissRes2 =  lgm(swissRain, locations=30, formula=rain ~ elev + factor(land),
		covariates=list(elev=swissAltitude,land=swissLandType), 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE
)
swissRes2$summary
bob(swissRes2)



# simulated data (without a CRS)
# and all covariates are in 'data' object
myModel = c(intercept=0,variance=2^2,nugget=0.5^2, range=2.5,shape=2, 
		cov1=0.2, cov2=-0.5)
covariates = brick(
		xmn=0,ymn=0,xmx=10,ymx=10,
		ncols=200,nrows=200,nl=2)
values(covariates)[,1] = rep(seq(0,1,len=nrow(covariates)), ncol(covariates))
values(covariates)[,2] = rep(seq(0,1,len=nrow(covariates)), 
		rep(nrow(covariates), ncol(covariates)))
names(covariates) = c("cov1","cov2")

Npoints = 40
myPoints = SpatialPoints(cbind(runif(Npoints,0,10), runif(Npoints,0,10)))	
# check for points too close together
thedist = spDists(myPoints)
thedist[lower.tri(thedist,diag=TRUE)]=NA
thedist = apply(thedist<0.2,2, any,na.rm=TRUE)
myPoints = myPoints[!thedist]


myPoints = SpatialPointsDataFrame(myPoints, 
		data=as.data.frame(extract(covariates, myPoints)))
library("RandomFields")
myPoints$U = geostatsp::RFsimulate(myModel,myPoints)$variable1 
myPoints$y= myModel["intercept"] +
		as.matrix(myPoints@data[,names(covariates)]) %*% 
		myModel[names(covariates)] +
		myPoints$U+
		rnorm(length(myPoints), 0, sqrt(myModel["nugget"]))

fitLikfit = likfitLgm(myPoints, trend=y~cov1+cov2, 
		param=c(range=1,nugget=0,shape=1)) 


Srange = exp(seq(log(0.75), log(6), len=20))

Slik = NULL
SlikWithN=NULL
Snugget=NULL
for(D in Srange) {
	Slik = c(Slik,
			loglikLgm(param=c(range=D,nugget=0,shape=1),
					data=myPoints, trend=y~cov1+cov2))
	temp = likfitLgm(paramToEstimate = "nugget",
			param=c(range=D,nugget=5,shape=1),
			data=myPoints, trend=y~cov1+cov2)
	SlikWithN = c(SlikWithN,
			temp$opt$value
			)	
	Snugget = c(Snugget, temp$opt$par["nugget"])		
}

pdf("profileL.pdf",height=4,width=12)
par(mfcol=c(1,3))
plot(Srange, Slik, ylab="-log lik", main="nugget=0")
plot(Srange, SlikWithN, ylab="-log lik", main="estimate nugget")
plot(Srange, Snugget, ylab="optimal nugget", main="the nugget")
dev.off()


# run lgm without providing covariates
fitMLE =  lgm(myPoints, locations=10, formula=y~ cov1 + cov2, 
		shape=1, fixShape=TRUE)


c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])
bob(fitMLE)

# now give covariates as raster brick
fitMLE =  lgm(myPoints, locations=10, formula=y~ cov1 + cov2, 
		covariates=covariates,
		shape=1, fixShape=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])
bob(fitMLE)
# now give covariates as list
fitMLE =  lgm(myPoints, locations=10, formula=y~ cov1 + cov2, 
		covariates=list(cov1=covariates[["cov1"]],
				cov2 = covariates[["cov2"]]),
		shape=1, fixShape=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])
bob(fitMLE)

# not remove covariates from data
myPoints = SpatialPointsDataFrame(SpatialPoints(myPoints),
		data=myPoints@data[,"y",drop=FALSE])

# now give covariates as raster brick
fitMLE =  lgm(myPoints, locations=10, formula=y~ cov1 + cov2, 
		covariates=covariates,
		shape=1, fixShape=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])
bob(fitMLE)

