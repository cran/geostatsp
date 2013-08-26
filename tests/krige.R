# land type, categorical variable
library("geostatsp")
data("swissRain")
swissRain$lograin = log(swissRain$rain)
swissRain$elevation = extract(swissAltitude, swissRain)
swissAltitude[1:50,1:50] = NA

swissRaster = raster(extent(swissBorder), ncols=20, nrows=20, 
		crs=swissRain@proj4string)	


swissRain$land = raster::extract(swissLandType, swissRain)
# get rid of land types with few observations
landTable = table(swissRain$land)
landTable = as.numeric(names(landTable)[landTable > 5])
swissRain2 = swissRain [swissRain$land %in% landTable, ]

swissFit3 = likfitLgm(data=swissRain2, 
		trend=lograin~ elevation + factor(land),
		param=c(range=46500, nugget=0.05,rough=1,  
				aniso.angle.degrees=35, aniso.ratio=12),
		paramToEstimate = c("range","nugget", 
				"aniso.angle.degrees", "aniso.ratio"),
		parscale = c(range=5000,nugget=0.01, 
				aniso.ratio=1,aniso.angle.degrees=5)
)

swissKrige3 = krige(data=swissRain2, trend = swissFit3$trend,
		param=swissFit3$param, 
		covariates = list(elevation = swissAltitude,land=swissLandType),
		locations = swissRaster, expPred=TRUE)

plot(swissKrige3[["predict"]])	
plot(swissBorder, add=TRUE)


# now change land to a factor
landTypes = swissLandType@data@attributes[[1]]
# remove land types with no observations on them
landTypes=landTypes[landTypes[,1] %in% unique(swissRain2$land),]

swissRain2$landFac = factor(swissRain2$land, 
		levels=landTypes[,1],
		labels=landTypes[,2])

swissFit4 = likfitLgm(data=swissRain2, 
		trend=lograin~ elevation + landFac,
		param=c(range=46500, nugget=0.05,rough=1,  
				aniso.angle.degrees=35, aniso.ratio=12),
		paramToEstimate = c("range","nugget", 
				"aniso.angle.degrees", "aniso.ratio"),
		parscale = c(range=5000,nugget=0.01, 
				aniso.ratio=1,aniso.angle.degrees=5)
)
swissKrige4 = krige(data=swissRain2, trend = swissFit4$trend,
		param=swissFit4$param, 
		covariates = list(elevation = swissAltitude,landFac=swissLandType),
		locations = swissRaster,expPred=TRUE )
plot(swissKrige4[["predict"]])	
plot(swissBorder, add=TRUE)

swissRain2$landFac2 = as.character(swissRain2$landFac)

swissFit5= likfitLgm(data=swissRain2, 
		trend=lograin~ elevation + factor(landFac2),
		param=c(range=46500, nugget=0.05,rough=1,  
				aniso.angle.degrees=35, aniso.ratio=12),
		paramToEstimate = c("range","nugget", 
				"aniso.angle.degrees", "aniso.ratio"),
		parscale = c(range=5000,nugget=0.01, 
				aniso.ratio=1,aniso.angle.degrees=5)
)


swissKrige5 = krige(data=swissRain2, trend = swissFit5$trend,
		param=swissFit5$param, 
		covariates = list(elevation = swissAltitude,landFac2=swissLandType),
		locations = swissRaster,expPred=TRUE)

plot(swissKrige5[["predict"]])	
plot(swissBorder, add=TRUE)

