library("geostatsp")
data("swissRain")


# specify formula name of raster layer
swissFit = lgm(data=swissRain, formula=rain~ SRTM_1km,
		locations=80, covariates=swissAltitude,
		rough=1,  fixRough=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFit)
swissFit$param


# specify formula using name of list element

swissFitAgain = lgm(data=swissRain, formula=rain~ elev,
		locations=80, covariates=list(elev=swissAltitude),
		rough=1,  fixRough=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param


swissFitAgain = lgm(data=swissRain, formula="rain",
		locations=80, covariates=swissAltitude,
		rough=1,  fixRough=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param
swissFitAgain = lgm(data=swissRain, formula="rain",
		locations=80, covariates=list(elev=swissAltitude),
		rough=1,  fixRough=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE)	
names(swissFitAgain)
swissFitAgain$param

# simulated data (without a CRS)
# and all covariates are in 'data' object
myModel = c(mean=0,variance=1,nugget=0, range=0.5,rough=2)

cov1 = cov2 = raster(extent(c(0,10,0,10)), ncol=100,nrow=100)
values(cov1) = rep(seq(0,1,len=nrow(cov1)), ncol(cov1))
values(cov2) = rep(seq(0,1,len=nrow(cov1)), rep(nrow(cov1), ncol(cov1)))
myU = GaussRF(cov1,param=myModel)
myLambda = 0.5 + 0.2*cov1 - 0.5 * cov2 + myU 

myPoints = SpatialPoints(cbind(runif(40,0,10), runif(40,0,10)))	
myObs = SpatialPointsDataFrame(myPoints,
		data=data.frame(y= extract(myLambda, myPoints) + 
						rnorm(length(myPoints), mean=0, sd=0.5) 
		)	) 

myCov = list(cov1 = cov1, cov2=cov2)

fitMLE =  lgm(myObs, locations=10, formula=y~ cov1 + cov2, 
		covariates = myCov,
		rough=1, fixRough=TRUE)
fitMLE$summary["range","estimate"]



