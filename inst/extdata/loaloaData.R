
	
#####################
# disease status
######################

	loaloa = read.table(
			"http://www.leg.ufpr.br/lib/exe/fetch.php/pessoais:paulojus:mbgbook:datasets:loaloa.txt",
			skip=7)
	names(loaloa) = c("long","lat", "N","y","elevation","ndvi","sdndvi")
	library(sp)
	library(rgdal)
	
# create projection without epsg code so rgdal doesn't need to be loaded
	theproj = CRS("+init=epsg:3064")
	
	loaloaLL = SpatialPointsDataFrame(loaloa[,c("long","lat")], 
			data=loaloa[,c("N","y")], 
			proj4string = CRS("+init=epsg:4326"))
	loaloa=spTransform(loaloaLL, theproj)
	loaloa$villageID = seq(1, length(loaloa))
	

	
	library('geostatsp')
	library('mapmisc')
	data('loaloa')
	loaM = spTransform(loaloa, crsModis)

	loaExtentM = extend(extent(loaM), 40*1000)
	loaExtent = extend(extent(loaloa), 20*1000)
	
	loaTiles = factorValues(modisRaster,
      values(raster::crop(
	      			modisRaster, extend(extent(loaM), 100*1000), snap='out'
      		)))[,'tile']
	loaTilesCollapse = paste(loaTiles, collapse='|')
	dataDir = "/store/patrick/spatialData/"
	
	#################
# land cover data
	#################
	
	myProduct = "MCD12Q1"
	
	modisUrl = 'ftp://ladsweb.nascom.nasa.gov/allData/5/MCD12Q1/2002/001/'
	
	
	Sfiles = paste(modisUrl,
  		grep(
					loaTiles,
    			unlist(strsplit(RCurl::getURL(
	    								modisUrl,ftp.use.epsv=TRUE,
	    								dirlistonly = TRUE), '\\n')
	  			),
  				value=TRUE), 
			sep='')
	
	
	Slocal = Pmisc::downloadIfOld(
			Sfiles,
			file.path(dataDir, basename(Sfiles))
	)	
	
	rList = mapply(
			function(x) stack(gdalUtils::get_subdatasets(x)),
			Slocal
	)
	
	landOrig = merge(rList[[1]][[1]], rList[[2]][[1]])
	
	landCrop = raster::crop(
			landOrig, loaExtentM
	)
	
	ltLoaExt = projectRaster(landCrop, crs=crs(loaloa), method='ngb')
	
	ltLoa = raster::crop(
			ltLoaExt, loaExtent
	)
	
	ltLoa = ratify(ltLoa)
	
	
	landTable = XML::readHTMLTable(RCurl::getURL(
					'https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd12q1'
			), header=TRUE, stringsAsFactors = FALSE)
	
	landTable = landTable[[2]][,1:2]
	
	
	landTableM = colourScale(
			x=ltLoa, labels=landTable, style='unique', breaks=10,
			col='Set3', exclude=0
	)
	
	ltLoa@legend@colortable = landTableM$colortable
	levels(ltLoa)[[1]] = landTableM$levels
	
	plot(ltLoa)
	legendBreaks("bottomleft", ltLoa, width=20, cex=0.7)
	points(loaloa)

#############3
# EVI
##############

	
# EVI, get 12 month average for 2002
	myProduct = "MOD13Q1"
	
	modisUrl = 'ftp://ladsweb.nascom.nasa.gov/allData/5/MOD13Q1/2002/'
	
		
	urlWithDate = paste(modisUrl, unlist(strsplit(RCurl::getURL(
	    				modisUrl,ftp.use.epsv=TRUE,
	    				dirlistonly = TRUE), '\\n')
	), '/', sep='')
	

	allDates = strsplit(RCurl::getURL(
	    	urlWithDate, ftp.use.epsv=TRUE,
	    	dirlistonly = TRUE), '\\n')

	allDates = paste(
			rep(names(allDates), unlist(lapply(allDates, length))),
			'/',
			unlist(allDates),
			sep=''
			)
	
	Sfiles= grep(loaTilesCollapse, allDates, value=TRUE)		

	Slocal = Pmisc::downloadIfOld(
			Sfiles,
			file.path(dataDir, basename(Sfiles))
	)	
	
	Sevi = mapply(
			function(x)
				grep('EVI$', gdalUtils::get_subdatasets(x), value=TRUE),
			x=Slocal)
	
	
	rList = mapply(
			function(x) {
				theFile = file.path(dataDir, paste(x, '.grd', sep=''))
				brick(stack(grep(x, Sevi, value=TRUE)), 
						filename = theFile, overwrite = file.exists(theFile))
				
			},
			x=loaTiles)

	mFile = file.path(dataDir, 'eviMerge.grd')
	eviModisTime = merge(
			crop(rList[[1]], loaExtentM),
			crop(rList[[2]], loaExtentM),
			file=mFile, overwrite=file.exists(mFile)
	)
	
	cFile = file.path(dataDir, 'eviClamp.grd')
	eviModisC= clamp(eviModisTime, 5*10^6, Inf, 
			useValues=FALSE,
			file=cFile, overwrite=file.exists(cFile))
	
	eviModisAgg = aggregate(eviModisC, fact=4, fun=mean)
	eviModis = calc(eviModisAgg, mean, na.rm=TRUE)
	
	eviExt = projectRaster(eviModis, crs=projection(loaloa))
	
	eviLoa = raster::crop(
			eviExt, loaExtent
	)
	names(eviLoa) = 'evi'
	
	eviCol = colourScale(eviLoa, dec=-7,
			col='Greens', breaks=12, style='equal')
	
	map.new(loaloa, buffer=30*1000)
	plot(eviLoa,add=TRUE, legend=FALSE,col=eviCol$col,breaks=eviCol$breaks)
	legendBreaks('bottomleft', eviCol)
	points(loaloa)


	
	############
# elevation
	#############
	myProduct = "SRTM"
	getProduct(myProduct)
	getCollection(myProduct)
# doesn't seem to be available!
# do it the hard way...
	
# download six tiles
	
	theTiles=c("3811","3911","3812","3912", "4011", "4012")
	if(FALSE) {
		for(D in theTiles) {
			D2 = paste(substr(D,1,2), "_", substr(D,3,4),sep="")
			print(D2)
			
			download.file(paste(
							"ftp://srtm.csi.cgiar.org/SRTM_V41/SRTM_Data_GeoTiff/srtm_",
							D2, ".zip",sep=""),
					paste(dataDir, "srtm_", D2, ".zip",sep="")
			)
		}
	}
	
# extent to crop to
	forCrop = extend(projectExtent(
					loaloa,crs="+init=epsg:4326")@extent,
			0.5)
	
	for(D in theTiles) {
		print(D)
		D2 = paste(substr(D,1,2), "_", substr(D,3,4),sep="")
		print(D2)
		
		unzip(paste(dataDir,"srtm_", D2, ".zip",sep=""), exdir=dataDir)
		elevUnc = raster(paste(dataDir, "srtm_", D2, ".tif",sep=""))
		
		elevC = crop(elevUnc, forCrop)
		
		if(D == theTiles[1]) {
			result = elevC
		} else {
			result = merge(result, elevC)
		}
	}
	
	plot(result)
	
	elevAgg=raster::aggregate(result, fact=30)
	
	
	elevationLoa = projectRaster(elevAgg, crs=CRS(proj4string(loaloa)))
	
  
  
	#################
# temperature
	###################
	
	myProduct = "MOD11C3"
	thehdf=getHdf(product=myProduct,begin="2002-01-01",end="2002-12-31", 
			extent=extent(spTransform(loaloa, CRS("+init=epsg:4326"))))
	
# the HDR files have multiple layers
# find the layer with EVI data
	layerNames = getSds(thehdf[[1]][1])$SDSnames
	theLayer = grep("^LST_Day", layerNames)
	theString = rep(0, length(layerNames))
	theString[theLayer] = 1
	theString = paste(theString, collapse="")
	
# convert downloaded HDR files to nice tiff files
# on a 2km grid, with the same projection as loaloa
	runGdal(product=myProduct,begin="2002-01-01",end="2002-12-31",
			SDSstring = theString, job="loa"
	) 
	
# find names of all the tiff files 
	thenames = preStack(
			path = paste(options()$MODIS_outDirPath, "loa",sep=""),
			pattern=myProduct)
# create a raster stack of EVI for each day
	tempLoa= stack(thenames)
# crop out the study region (20km buffer around the loaloa data)
	tempLoa = crop(tempLoa, 
			extend(extent(projectExtent(loaloa,crs="+init=epsg:4326")),0.5)
	)
	tempLoa = projectRaster(tempLoa, crs=proj4string(loaloa))
	
# compute the yearly average
	tempAvg = raster::overlay(tempLoa, fun=function(x) {
				mean(x, na.rm=T)
			})
	
# modis says the temperatures are scaled by 0.02
# to get degrees celcius do
	tempLoa = tempAvg*0.02-273.15
# plot the data
	plot(tempLoa)
	points(loaloa)
	range(extract(tempLoa, loaloa))
  
	
	
# save data
	save(loaloa, eviLoa, ltLoa, tempLoa, elevationLoa,
			file="~/research/diseasemapping/pkg/geostatsp/data/loaloa.RData", 
			compress="xz")
	
	

