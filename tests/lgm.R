library("geostatsp")
data("swissRain")
swissRain = unwrap(swissRain)
swissAltitude= unwrap(swissAltitude)
swissLandType = unwrap(swissLandType)

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
swissFit = lgm(data=swissRain[1:60,], formula=rain~ CHE_alt,
		grid=80, covariates=swissAltitude,
		shape=2,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFit)
swissFit$summary[,1:4]
bob(swissFit)


# specify formula using name of list element

swissFitAgain = lgm(
    data=swissRain[1:80,], 
    formula=rain~ elev+land,
		grid=80, covariates=list(elev=swissAltitude,land=swissLandType),
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE,
		aniso=TRUE, verbose=FALSE)	
names(swissFitAgain)
swissFitAgain$param
swissFitAgain$summary[,c('estimate','stdErr','Estimated')]
bob(swissFitAgain)

swissFitAgain = lgm(
	data=swissRain[1:60,], 
	formula="rain",
		grid=80, covariates=swissAltitude,
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param
bob(swissFitAgain)


  

swissFitAgain = lgm(data=swissRain, formula="rain",
		grid=80, covariates=list(elev=swissAltitude,land=swissLandType),
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param
bob(swissFitAgain)


# land type, factor covariate
swissRes2 =  lgm(rain ~ elev + factor(land), swissRain, 
		grid=30, 
		covariates=list(elev=swissAltitude,land=swissLandType), 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE
)
swissRes2$summary
bob(swissRes2)


  
# simulated data (without a CRS)
# and all covariates are in 'data' object
myModel = c(intercept=0,variance=2^2,nugget=0.5^2, range=4.5,shape=2, 
		cov1=0.2, cov2=-0.5)
covariates = rast(
		xmin=0,ymin=0,xmax=10,ymax=10,
		ncols=200,nrows=200,nlyrs=2, crs='+proj=merc')
values(covariates) = c(
		rep(seq(0,1,len=nrow(covariates)), ncol(covariates)),
		rep(seq(0,1,len=nrow(covariates)), 
		rep(nrow(covariates), ncol(covariates))))
names(covariates) = c("cov1","cov2")

Npoints = 30
set.seed(0)
myPoints = vect(cbind(runif(Npoints,0,10), 
        seq(0,10, len=Npoints)), crs = '+proj=merc')

# don't test using the randomFields package, it's currently broken on some R builds
options(useRandomFields = FALSE)

myPoints = RFsimulate(myModel,myPoints)

values(myPoints) = cbind(
   values(myPoints), 
  as.data.frame(extract(covariates, myPoints))
)

myPoints$y= myModel["intercept"] +
		as.matrix(values(myPoints)[,names(covariates)]) %*% 
		myModel[names(covariates)] +
		myPoints$sim+
		rnorm(length(myPoints), 0, sqrt(myModel["nugget"]))

fitLikfit = likfitLgm(y~cov1+cov2, myPoints,  
		param=c(range=1,nugget=0,shape=1)) 




# run lgm without providing covariates
fitMLE =  lgm(
    formula=y~ cov1+cov2, 
    data=myPoints,
    grid=10, covariates=list(), 
		shape=1, fixShape=TRUE)

c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])
bob(fitMLE)


# now give covariates as raster brick
fitMLE =  lgm( y~ cov1 + cov2, myPoints, grid=10,  
		covariates=covariates,
		shape=1, fixShape=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])
bob(fitMLE)
# now give covariates as list
fitMLE =  lgm(y~ cov1 + cov2, myPoints, grid=10,   
		covariates=list(cov1=covariates[["cov1"]],
				cov2 = covariates[["cov2"]]),
		shape=1, fixShape=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])
bob(fitMLE)

# not remove covariates from data
myPoints = vect(crds(myPoints),
		atts=values(myPoints)[,"y",drop=FALSE])

# now give covariates as raster brick
fitMLE =  lgm(y~ cov1 + cov2,  myPoints, grid=10,  
		covariates=covariates,
		shape=1, fixShape=TRUE)
c(fitMLE$summary["range","estimate"], fitLikfit$summary["range","estimate"])
bob(fitMLE)


