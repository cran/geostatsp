\name{swissRain}
\alias{swissRain}
\alias{swissAltitude}
\alias{swissBorder}
\alias{swissLandType}

\docType{data}
\title{
Swiss rainfall data
}
\description{
Data from the SIC-97 project: Spatial Interpolation Comparison.
}
\usage{data("swissRain")}
\format{
\code{swissRain} is a SpatialPolygonsDataFrame 
100 daily rainfall
measurements made in Switzerland on the 8th of May 1986.
\code{swissAltitude} is a raster of elevation data, and \code{swissLandType} is a raster
of land cover types.
}
\source{
\url{https://wiki.52north.org/AI_GEOSTATS/AI_GEOSTATSData} and  \url{https://srtm.csi.cgiar.org} and
\url{https://lpdaac.usgs.gov/product_search/?collections=Combined+MODIS&collections=Terra+MODIS&collections=Aqua+MODIS&view=list}
}

\examples{
data("swissRain")
plot(swissAltitude, main="elevation")
points(swissRain)
plot(swissBorder, add=TRUE)


# land type, a categorical variable
commonValues  = sort(table(values(swissLandType)),decreasing=TRUE)[1:5]
commonValues=commonValues[!names(commonValues)==0]

thelevels = levels(swissLandType)[[1]]$ID
thebreaks = c(-0.5, 0.5+thelevels)
thecol = rep(NA, length(thelevels))
names(thecol) = as.character(thelevels)

thecol[names(commonValues)] = rainbow(length(commonValues))

plot(swissLandType, breaks=thebreaks, col=thecol,legend=FALSE,
	main="land type")
points(swissRain)
plot(swissBorder, add=TRUE)


legend("topleft",fill=thecol[names(commonValues)],
		legend=levels(swissLandType)[[1]][
						match(as.integer(names(commonValues)),
								levels(swissLandType)[[1]]$ID),
						"Category"],
				bty='n'
				)

# code to assemble the dataset
\dontrun{
	dataDir = tempdir()

	download.file(
"https://wiki.52north.org/pub/AI_GEOSTATS/AI_GEOSTATSData/sic97data_01.zip",
		destfile=paste(dataDir, "swiss.zip",sep=""))
	swissFile = unzip(paste(dataDir, 'swiss.zip',sep=""), exdir=dataDir)
	swissRain = read.table(grep("sic_obs.dat",swissFile, value=TRUE), 
		sep=',', 
		col.names=c('ID','x','y','rain'))
	# the following seems to make the coordinates line up with epsg:2056
	swissRain$x = swissRain$x - 17791.29 + 2672591 
	swissRain$y = swissRain$y - 13224.66 + 1200225
	# the readme file says rain is in tenths of mm 
	swissRain$rain= swissRain$rain / 10  
	library(sp)
	library(rgdal)
	# create projection without epsg code so rgdal doesn't need to be loaded
	#theproj = CRSargs(CRS("+init=epsg:2056"))
	theproj = CRS(paste("+proj=somerc +lat_0=46.9524055555556",
		"+lon_0=7.43958333333333 +k_0=1 +x_0=2600000 +y_0=1200000"))

	#theproj = gsub("\\+init=epsg:[[:digit:]]+ ", "", theproj)
	#theproj = CRS(theproj)

	swissRain = SpatialPointsDataFrame(swissRain[,c('x','y')], data=swissRain[,c('ID','rain')], 
		proj4string=theproj)



	#######
	# Swiss Border
	#######


	library('raster')
	swissBorder = getData('GADM', country='CHE', level=0)
	isChar = which(unlist(lapply(swissBorder@data, is.character)))
	isUtf = which(
		unlist(lapply(swissBorder@data[,isChar], 
			Encoding)) == 'UTF-8')
	swissBorder = swissBorder[,
		!match(names(swissBorder), names(isUtf), nomatch=0)
		]
	library(rgdal)
	swissBorder = spTransform(swissBorder, CRS(proj4string(swissRain)))

############
# land type
############
# see loaloa's help file for installation of the MODIS package
library(MODIS)
MODISoptions(gdalPath="/usr/bin/", 
	localArcPath=dataDir, outDirPath=dataDir)
options()[grep("MODIS", names(options()), value=TRUE)]

myProduct = "MCD12Q1"
getProduct(myProduct)

thehdf=getHdf(product=myProduct,
		begin="2002-01-01",end="2002-01-02",
		tileH=18, tileV=4)
#		extent=extent(spTransform(swissBorder, CRS("+init=epsg:4326"))))

layerNames = getSds(thehdf[[1]][1])$SDSnames
ltLayer = grep("Type_1$", layerNames)
theString = rep(0, length(layerNames))
theString[ltLayer] = 1
theString = paste(theString, collapse="")

runGdal(product=myProduct,
		begin="2002-01-01",end="2002-01-02",
		outProj = proj4string(swissRain),
		pixelSize=2000, job="loa",
		SDSstring = theString,
				tileH=17:18, tileV=3:4)

#		extent=extent(spTransform(swissBorder, CRS("+init=epsg:4326"))))

# find file name 
thenames = preStack(
		path = paste(options()$MODIS_outDirPath, "loa",sep=""),
		pattern=myProduct)
swissLandType = raster(thenames)
swissLandType = crop(swissLandType, extend(extent(swissBorder),20000))


swissLandType = as.factor(swissLandType)

# labels of land types
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


levels(swissLandType)[[1]]$Category = 
		labelVec[as.character(levels(swissLandType)[[1]]$ID)]

levels(swissLandType)[[1]]$col = NA

theForests = grep("forest", levels(swissLandType)[[1]]$Category, 
	ignore.case=TRUE)
	
library(RColorBrewer)	
levels(swissLandType)[[1]][theForests,"col"] = 
	brewer.pal(length(theForests)+1, "Greens")[-1]
	
levels(swissLandType)[[1]][
	grep("snow", levels(swissLandType)[[1]]$Category,ignore.case=TRUE),
	"col"] = "#FFFFFF"

levels(swissLandType)[[1]][
	grep("water", levels(swissLandType)[[1]]$Category,ignore.case=TRUE),
	"col"] = "#0000FF"

levels(swissLandType)[[1]][
	grep("grass", levels(swissLandType)[[1]]$Category,ignore.case=TRUE),
	"col"] = "#CCBB00"

	
stillNA = is.na(levels(swissLandType)[[1]]$col)
levels(swissLandType)[[1]][stillNA, "col"] =
	brewer.pal(sum(stillNA), "Set3")
	
	

swissLandType@legend@colortable = levels(swissLandType)[[1]]$col

levels(swissLandType)[[1]]$n = table(values(swissLandType))

plot(swissLandType)
mostCommon = levels(swissLandType)[[1]]$n >= 700
legend("topright", 	
	fill=levels(swissLandType)[[1]][mostCommon,"col"],
	legend = substr(
		levels(swissLandType)[[1]][mostCommon,"Category"],
		1, 12)
	)
	
table(extract(swissLandType, swissRain), exclude=NULL)


#### 
# SwissAltitude
###
	library(raster)
	download.file('http://biogeo.ucdavis.edu/data/diva/alt/CHE_alt.zip',
		destfile=paste(dataDir,'CHE_alt.zip',sep=""))
	swissAfile = unzip(paste(dataDir, 'CHE_alt.zip',sep=""), exdir=dataDir)
	swissAltitude = raster(grep("CHE_alt.gri", swissAfile, value=TRUE))
		swissAltitude = projectRaster(swissAltitude, 
			crs=CRS(proj4string(swissRain)))
		swissAltitude = aggregate(swissAltitude,fact=2)
	 
 
save(swissRain, swissAltitude, swissBorder, swissLandType,
		file="~/research/diseasemapping/pkg/geostatsp/data/swissRain.RData",
		compress="xz")
}

}


\keyword{datasets}
