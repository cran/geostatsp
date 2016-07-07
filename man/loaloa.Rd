\name{loaloa}
\alias{loaloa}
\alias{elevationLoa}
\alias{eviLoa}
\alias{ltLoa}
\alias{tempLoa}
\docType{data}
\title{
Loaloa prevalence data from 197 village surveys
}
\description{
Location and prevalence data from villages, elevation an vegetation index for the study region.
}
\usage{data("loaloa")}
\format{
\code{loaloa} is a SpatialPolygonsDataFrame of the data, with columns \code{N} being the number 
of individuals tested and \code{y} being the number of positives.
\code{elevationLoa} is a raster of elevation data.
\code{eviLoa} is a raster of vegetation index for a specific date. \code{ltLoa} is land type.
\code{ltLoa} is a raster of land types. 1     2     5     6     7     8     9    10    11    12    13    14    15
\code{tempLoa} is a raster of average temperature in degrees C.
}

\source{
\url{http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets} for the loaloa data,
\url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1} for the EVI data,  
\url{https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd12q2} for land type and  
\url{http://srtm.csi.cgiar.org} for the elevation data.
}

\examples{


data("loaloa")
plot(loaloa, main="loaloa villages")

# elevation
plot(elevationLoa, col=terrain.colors(100), main="elevation")
points(loaloa)

# vegetation index
plot(eviLoa, main="evi")
points(loaloa)

plot(tempLoa, main="temperature")
points(loaloa)



# land type, a categorical variable
plot(ltLoa)
text(768000, 800000,"land type")
commonValues  = order(levels(ltLoa)[[1]]$n,decreasing=TRUE)[1:7]
# get rid of water
commonValues = commonValues[commonValues != 1]

legend("bottomleft",
	fill=levels(ltLoa)[[1]][commonValues,"col"],
		legend=levels(ltLoa)[[1]][commonValues,"Category"],
				bty='n'
)

points(loaloa)

# the following is the code for creating these datasets
\dontrun{
dataDir = "/store/patrick/spatialData/"

loaloa = read.table(
"http://www.leg.ufpr.br/lib/exe/fetch.php/pessoais:paulojus:mbgbook:datasets:loaloa.txt",
skip=7)
names(loaloa) = c("long","lat", "N","y","elevation","ndvi","sdndvi")
library(sp)
library(rgdal)

# create projection without epsg code so rgdal doesn't need to be loaded
theproj = CRSargs(CRS("+init=epsg:3064"))
theproj = gsub("\\+init=epsg:[[:digit:]]+ ", "", theproj)
theproj = CRS(theproj)

loaloaLL = SpatialPointsDataFrame(loaloa[,c("long","lat")], 
		data=loaloa[,c("N","y")], 
		proj4string = CRS("+init=epsg:4326"))
loaloa=spTransform(loaloaLL, theproj)
loaloa$villageID = seq(1, length(loaloa))



if(FALSE) {
	install.packages("RCurl")
	install.packages("rgeos")
	install.packages('mapdata')
	install.packages("MODIS", repos="http://r-forge.r-project.org") 
}
	
# load the package and set options

library(MODIS)
MODISoptions(gdalPath="/usr/bin/", localArcPath=dataDir)
options()[grep("MODIS", names(options()), value=T)]

# EVI, get 12 month average for 2002
myProduct = "MOD13Q1"
getProduct(myProduct)

thehdf=getHdf(product=myProduct,begin="2002-01-01",end="2002-12-31",
		extent=extent(spTransform(loaloa, CRS("+init=epsg:4326"))))

# the HDR files have multiple layers
# find the layer with EVI data
layerNames = getSds(thehdf[[1]][1])$SDSnames
eviLayer = grep("EVI$", layerNames)
theString = rep(0, length(layerNames))
theString[eviLayer] = 1
theString = paste(theString, collapse="")

# convert downloaded HDR files to nice tiff files
# on a 2km grid, with the same projection as loaloa
runGdal(product=myProduct,begin="2002-01-01",end="2002-12-31",
		outProj = proj4string(loaloa),
		pixelSize=2000, job="loa",
		SDSstring = theString,
		extent=extent(spTransform(loaloa, CRS("+init=epsg:4326")))) 

# find names of all the tiff files 
thenames = preStack(
		path = paste(options()$MODIS_outDirPath, "loa",sep=""),
		pattern=myProduct)
# create a raster stack of EVI for each day
eviLoa= stack(thenames)
# crop out the study region (20km buffer around the loaloa data)
eviLoa = crop(eviLoa, 
		extend(extent(loaloa),20000))

# compute the yearly average
eviAvg = raster::overlay(eviLoa, fun=function(x) {
			mean(x, na.rm=T)
		})

# plot the data
plot(eviAvg)
points(loaloa)
range(extract(eviLoa, loaloa))
# save data as an R object
eviLoa = eviAvg

#################
# land cover data
#################
myProduct = "MCD12Q1"
getProduct(myProduct)

thehdf=getHdf(product=myProduct,
		begin="2001-01-01",end="2010-12-30",
		extent=extent(
				spTransform(loaloa, CRS("+init=epsg:4326"))))

layerNames = getSds(thehdf[[1]][1])$SDSnames
ltLayer = grep("Type_1$", layerNames)
theString = rep(0, length(layerNames))
theString[ltLayer] = 1
theString = paste(theString, collapse="")

runGdal(product=myProduct,
		begin="2002-01-01",end="2002-12-30",
		outProj = proj4string(loaloa),
		pixelSize=2000, job="loa",
		SDSstring = theString,
		extent=extent(spTransform(loaloa, CRS("+init=epsg:4326"))))

# find file name 
thenames = preStack(
		path = paste(options()$MODIS_outDirPath, "loa",sep=""),
		pattern=myProduct)
ltLoa = raster(thenames)
ltLoa = crop(ltLoa, extend(extent(loaloa),20000))

# convert to a raster of factors
ltLoa = as.factor(ltLoa)

# get land type labels
library(XML)
labels = readHTMLTable("http://nsidc.org/data/ease/ancillary.html")
labels = labels[[grep("Land Cover Classes", names(labels))]]
classCol = grep("Class Number", names(labels))
labels[,classCol] = as.integer(as.character(labels[,classCol]))

labels[  grep("Water", labels$Category), 
	classCol		
	] = 0
labelVec = as.character(labels$Category)
names(labelVec) = as.character(labels[,classCol])


levels(ltLoa)[[1]]$Category = 
		labelVec[as.character(levels(ltLoa)[[1]]$ID)]

# assign colours to most common land types
commonValues  = sort(table(values(ltLoa)),decreasing=TRUE)
ctable = levels(ltLoa)[[1]]
rownames(ctable) = as.character(ctable$ID)
ctable$n = NA
ctable[names(commonValues),"n"] = commonValues

rownames(ctable) = ctable$Category

ctable$col = NA
ctable["Snow and ice","col"] = "#FFFFFF"
ctable["Water bodies","col"] = NA
ctable["Urban and built-up","col"] = "#C31400"
ctable["Croplands","col"] = "#FFFF88"
ctable["Grasslands","col"] = "#CCDD11"
ctable["Permanent Wetlands","col"] ="#99BB99"
ctable["Cropland/natural vegetation mosaic","col"] = "#CC7711"
ctable["Barren or sparsely vegetated","col"] = "#111111"


forests = grep("forest",ctable$Category,value=T,ignore.case=T)
ctable[forests,"col"] = RColorBrewer::brewer.pal(length(forests), "Greens")

savannas = grep("savannas|shrublands",ctable$Category,value=T,ignore.case=T)
ctable[savannas,"col"] = RColorBrewer::brewer.pal(length(savannas), "Oranges")

# put labels back in the raster
levels(ltLoa) = list(ctable)

# assign colours to values
forctable = rep(NA, nrow(ctable))
forctable[ctable$ID+1] = ctable$col

ltLoa@legend@colortable = forctable


# plot land type
plot(ltLoa, legend=FALSE)
legend("bottomleft",
		fill=levels(ltLoa)[[1]]$col,
		legend=levels(ltLoa)[[1]][
						match(as.integer(names(commonValues)),
								levels(ltLoa)[[1]]$ID),
						"Category"]
	)

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
	file="~/workspace/diseasemapping/pkg/geostatsp/data/loaloa.RData", 
	compress="xz")


}

}


\keyword{datasets}
