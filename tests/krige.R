# land type, categorical variable

library("geostatsp")
data("swissRain")
swissRain = unwrap(swissRain)
swissAltitude = unwrap(swissAltitude)
swissLandType = unwrap(swissLandType)
swissBorder = unwrap(swissBorder)

swissRain$lograin = log(swissRain$rain)
swissRain$elevation = extract(swissAltitude, swissRain, ID=FALSE)
swissAltitude[1:50,1:50] = NA

swissRaster = rast(extent=ext(swissBorder), ncols=20, nrows=20, 
		crs=crs(swissRain))	


swissRain$land = extract(swissLandType, swissRain, ID=FALSE)
# get rid of land types with few observations
landTable = table(swissRain$land)
landTable = as.character(names(landTable)[landTable > 5])
swissRain2 = swissRain[as.character(swissRain$land) %in% landTable, ]
swissRain2$land = factor(swissRain2$land,
	levels = names(sort(table(swissRain2$land), decreasing=TRUE))[1:5])

swissFit3 = likfitLgm(
    data=swissRain2, 
		formula=lograin~ elevation + land,
		param=c(range=46500, nugget=0.05,shape=1,  
				anisoAngleDegrees=35, anisoRatio=12),
		paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"),
		parscale = c(range=5000,nugget=0.01, 
				anisoRatio=1,anisoAngleDegrees=5)
)

swissKrige3 = krigeLgm(
    data=swissRain2[1:60,], 
		formula = swissFit3$model$formula,
		param=swissFit3$param, 
		covariates = list(elevation = swissAltitude,land=swissLandType),
		grid = swissRaster, expPred=TRUE)

pdf("krige3.pdf")
plot(swissKrige3[["predict"]])	
plot(swissBorder, add=TRUE)
dev.off()


  

swissRain2$landFac2 = as.character(swissRain2$land)

swissFit5= likfitLgm(lograin~ elevation + landFac2,
		data=swissRain2, 
		param=c(range=46500, nugget=0.05,shape=1,  
				anisoAngleDegrees=35, anisoRatio=12),
		paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"),
		parscale = c(range=5000,nugget=0.01, 
				anisoRatio=1,anisoAngleDegrees=5)
)


swissKrige5 = krigeLgm(data=swissRain2[1:60,], 
		formula = swissFit5$model$formula,
		param=swissFit5$param, 
		covariates = list(elevation = swissAltitude,landFac2=swissLandType),
		grid = swissRaster,expPred=TRUE)
pdf("krige5.pdf")
plot(swissKrige5[["predict"]])	
plot(swissBorder, add=TRUE)
dev.off()






