library("geostatsp")
library("RandomFields")
model <- c(var=5, range=20,shape=0.5)

myraster = raster(nrows=100,ncols=100,xmn=0,ymn=0,xmx=10,ymx=10, 
		crs="+init=epsg:2081")

for(Dn in c(1,3)) {
	set.seed(0) 
	simu <- geostatsp::RFsimulate(model, x=myraster, n=Dn)
	set.seed(0) 
	simu2 <- geostatsp::RFsimulate(model, x=as(myraster,"SpatialPixels"), n=Dn)
	
	print(proj4string(simu))
	print(proj4string(simu2))
	
	par(mfrow=c(nlayers(simu),2))
	for(D in 1:nlayers(simu)) {
		plot(simu[[D]])
		plot(raster(simu2,layer=D))
	}
}

data("swissRain")
swissRain$sqrtrain = sqrt(swissRain$rain)

# estimate parameters


# isotropic
	swissRes =  lgm(swissRain, locations=20, formula="sqrtrain",
			covariates=swissAltitude,   
			shape=1, fixShape=TRUE,
			aniso=FALSE, fixNugget=FALSE,
			nuggetInPrediction=FALSE
	)
	
	
 # anisotropic
swissRes =  lgm(swissRain, locations=20, formula="sqrtrain",
		covariates=swissAltitude,   
		shape=1, fixShape=TRUE,
		aniso=TRUE, fixNugget=FALSE,
		nuggetInPrediction=FALSE
)
 


	# uncoinditional simulation
library("RandomFields")
RFoptions(printlevel=0)
swissSim = geostatsp::RFsimulate(
		model=swissRes$param,
		x=swissRes$predict,
		n=3
)


# simulate from the random effect conditional on
#   the observed data
swissSim = geostatsp::RFsimulate(
		model=swissRes$param,
		data=swissRes$resid,
		x=swissRes$predict,
		err.model=swissRes$param["nugget"],
		n=3
)

# plot the simulated random effect
plot(swissSim[[1]])
plot(swissBorder, add=TRUE)


# now with multiple parameter sets 
swissSim = geostatsp::RFsimulate(model=
				rbind(
						swissRes$param,
						swissRes$param*0.99),
		data=swissRes$resid,
		x=swissRes$predict,
		err.model=c(1, 0.99)*swissRes$param["nugget"],
		n=3
)
# plot the simulated random effect
plot(swissSim[[1]])
plot(swissBorder, add=TRUE)

# and multiple simulations
# now with multiple datasets 
swissSim = geostatsp::RFsimulate(model=
				rbind(
						swissRes$param,	
						0.99*swissRes$param,
						1.01*swissRes$param),
		data=swissRes$resid[,c(1,1,1)],
		err.model=c(1, 0.99, 1.01)*swissRes$param["nugget"],
		x=swissRes$predict,
		n=3
)
# plot the simulated random effect
plot(swissSim[[1]])
plot(swissBorder, add=TRUE)

# create a small raster of elevation data
swissAltSmall = aggregate(swissAltitude,fact=5)

# calculate the fixed effects portion of the rainfall process
rainMean = swissRes$param["(Intercept)"] +
		swissRes$param[ "CHE_alt" ] * swissAltSmall

# define a function to identify the location of maximum rainfall	
maxRainLocation = function(x, ...) {
	rain =  (rainMean + x)^2
	xyFromCell(rain, which.max(rain))
}


swissLocation = raster::cellStats(swissSim,   maxRainLocation)
plot(swissRes$predict[["predict"]])
plot(swissBorder, add=TRUE)
points(swissLocation)
