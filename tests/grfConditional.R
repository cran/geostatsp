library("geostatsp")
data("swissRain")
swissRain$elevation = extract(swissAltitude, swissRain)
swissRain$sqrtrain = sqrt(swissRain$rain)

# estimate parameters
swissFit =  likfitLgm(swissRain, trend=sqrtrain ~ elevation, 
		param=c(range=51700, nugget=0.11,rough=1,  
				aniso.angle.degrees=37, aniso.ratio=7.6),
		paramToEstimate = c("range","nugget", 
				"aniso.angle.degrees", "aniso.ratio"))

# simulate from the random effect conditional on
#   the observed data
swissSim = grfConditional(data=swissRain, ycol=swissFit$resid,
		param=swissFit$param, locations=20, 
		Nsim=1)

# plot the simulated random effect
pdf("grfConditional1.pdf")
plot(swissSim)
plot(swissBorder, add=TRUE)

# create a small raster of elevation data
swissAltSmall = aggregate(swissAltitude,fact=5)
# calculate the fixed effects portion of the rainfall process
rainMean = swissFit$param["(Intercept)"] +
			swissFit$param["elevation"] * swissAltSmall

# define a function to identify the location of maximum rainfall	
maxRainLocation = function(x) {
	rain =  (rainMean + x)^2
	xyFromCell(rain, which.max(rain))
}

# get a conditional sample of three locations of maximum rainfall
swissRain$resid = swissFit$resid
swissLocation = grfConditional(data=swissRain, 
		ycol="resid",
		param=swissFit$param, locations=swissAltSmall, 
		Nsim=3, fun = maxRainLocation)

# convert result from a list to a matrix
swissLocation = matrix(unlist(swissLocation), ncol=2,byrow=TRUE)
# add the locations to the map
points(swissLocation, pch=1:(dim(swissLocation)[1]),col='red')
dev.off()
