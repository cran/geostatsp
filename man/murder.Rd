\name{murder}
\alias{murder}
\alias{torontoPdens}
\alias{torontoIncome}
\alias{torontoNight}
\alias{torontoBorder}
\docType{data}
\title{
Murder locations
 }
\description{
  Locations of murders in Toronto 1990-2014
}
\usage{data("murder")}
\format{
  \code{murder} is a  SpatialPoints object of murder locations.  \code{torontoPdens}, 
  \code{torontoIncome}, and \code{torontoNight} are rasters containing
  population density (per hectare), median household income, and ambient light 
  respectively.  \code{torontoBorder} is a SpatialPolygonsDataFrame of the boundary of
  the city of Toronto.

}

\source{
Murder data:\url{https://www.thestar.com/news/crime/torontohomicidemap.html}, 

Lights: \url{https://ngdc.noaa.gov/eog/viirs/download_ut_mos.html}

Boundary files: \url{https://www150.statcan.gc.ca/n1/en/catalogue/92-160-X}

Income: \url{https://www150.statcan.gc.ca/n1/en/catalogue/97-551-X2006007}
}

\examples{
data("murder")
plot(torontoBorder)
points(murder, col="#0000FF40", cex=0.5)

data("torontoPop")

# light
data("murder")
plot(torontoNight, main="Toronto ambient light")
plot(torontoBorder, add=TRUE)
points(murder, col="#0000FF40", cex=0.5)


# income
plot(torontoIncome, main="Toronto Income")
points(murder, col="#0000FF40", cex=0.5)
plot(torontoBorder, add=TRUE)

# population density
plot(torontoPdens, main="Toronto pop dens")
points(murder, col="#0000FF40", cex=0.5)
plot(torontoBorder, add=TRUE)


\dontrun{
#building the dataset
fpath <- system.file("extdata", "murder1990.csv", package="geostatsp")
murderList=list()
# Load in datafiles retrieved from 
# http://www.thestar.com/news/crime/torontohomicidemap.html 
# Each year's murders are in a separate file, with 
# 1997, for example, being '/inst/extdata/murder1997.csv'
# this file was obtained by selecting the year '1997' from the
# menu marked 'select year', then clicking 'Download' at the bottom right
# selecting 'data' and in the new window clicking 'Underlying',
# then 'show all columnns' and 'dowload all rows as text file' 
for(Dyear in 1990:2014){
	Dfile = gsub("1990", Dyear, fpath)
	murderList[[as.character(Dyear)]] = read.table(Dfile, sep=",", header=TRUE,
		comment.char="", as.is=TRUE, quote="\"")
}
murderFull = do.call(rbind, murderList)

murder = murderFull[,c(
	"age.of.victim","sex.of.victim",
	"Homicide.type.groupings","method",
	 "URL","Details..if.available.")]
names(murder) = tolower(gsub("\\..+$", "", names(murder)))
names(murder) = gsub("^homicide$", "type", names(murder))
murder$method = gsub(",", "", murder$method)
murder$sex = factor(murder$sex)
murder$type = factor(murder$type)
murder$date = as.Date(
	as.character(murderFull$Date), 
	format="%B %d, %Y", tz='est')
murder$year = as.integer(format(murder$date, format="%Y"))
murder$date = as.Date(murder$date)

mcoord = as.matrix(
	murderFull[,c("Addresses._longitude",
		"Addresses._latitude")]
	)
	
murder=SpatialPointsDataFrame(
	coords=mcoord,
	data=murder, 
	proj4string=CRS("+init=epsg:4326")
	)
	
# get rid of children to keep this from being too morbid
murder = murder[which(murder$age >= 18),]

utm = CRS("+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
murder = spTransform(murder, utm)

# boundary of toronto
zipfileB = paste(tempfile(), ".zip", sep="")
download.file(
	paste('http://www12.statcan.gc.ca/census-recensement/2011/geo/',
	'bound-limit/files-fichiers/gcd_000b06a_e.zip',sep=''),
	zipfileB)
unzip(zipfileB, exdir=tempdir())
theshpB = grep("shp$",unzip(zipfileB,list=TRUE)$Name,  value=TRUE)
theshpB = gsub("\\.shp$", "", theshpB)
torontoBorder = rgdal::readOGR(tempdir(),theshpB,
	stringsAsFactors=FALSE)
torontoBorder = torontoBorder[
	grep("toronto", 
	as.character(torontoBorder$CDNAME), 
	ignore.case=TRUE),
]

 
torontoBorder =	spTransform(
	torontoBorder, CRS(proj4string(murder))
	)
library('mapmisc')
toMap = openmap(torontoBorder)
map.new(torontoBorder)
plot(toMap,add=TRUE)
plot(torontoBorder,add=TRUE)
points(murder,col='green')


save(murder, torontoBorder, 
	file="~/workspace/diseasemapping/pkg/geostatsp/data/murder.RData",
	compress='xz')

# census tract boundaries
zipfileCT = file.path(tempdir(),'toCt.zip')
download.file(
	paste('http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/',
	'files-fichiers/gct_000b06a_e.zip', sep=''),
		zipfileCT
		) 
		
unzip(zipfileCT, exdir=tempdir())
theshpCT = grep("shp$",unzip(zipfileCT,list=TRUE)$Name,  value=TRUE)
theshpCT = gsub("\\.shp$", "", theshpCT)

toCt2006 = readOGR(tempdir(),theshpCT)
toCt2006 = spTransform(toCt2006, CRS(proj4string(torontoBorder)))
aspoints = SpatialPoints(toCt2006)
projection(aspoints) = CRS(proj4string(torontoBorder))
inToronto = over(aspoints, torontoBorder)$CDNAME
toCt2006 =toCt2006[!is.na(inToronto),]
toCt2006 = spTransform(toCt2006, 
	CRS(proj4string(murder)))


# income
zipfileI = paste(tempfile(), ".zip", sep="")
download.file(
paste('http://www12.statcan.gc.ca/census-recensement/2006/dp-pd/tbt/',
'OpenDataDownload.cfm?PID=96273', sep=""),
zipfileI)

	unzip(zipfileI, exdir=tempdir())
	thesdmx = grep("^Generic",unzip(zipfileI,list=TRUE)$Name,  value=TRUE)
	library('rsdmx')
	toIncSdmx = readSDMX(file.path(tempdir(), thesdmx), isURL=FALSE)
	toInc = as.data.frame(toIncSdmx)
	toIncSub = toInc[grep("^535", toInc$GEO),]
	toIncSub = toIncSub[toIncSub$DIM0==3 & toIncSub$DIM1==1,]
	rownames(toIncSub) = toIncSub$GEO
toCt2006$income_median_household = toIncSub[
		gsub("\\.", "", as.character(toCt2006$CTUID)),
		'obsValue']

torontoIncome = rasterize(toCt2006, 
	squareRaster(torontoBorder, 150),
	field='income_median_household')



nightFile = paste(tempfile(), ".tif.gz", sep="")

# higher resolution at 
paste('http://mapserver.ngdc.noaa.gov/viirs_data/viirs_composite/'
'npp_20120418to20120426_20121011to20121023_sloff_15asec.'
'x2y1.75N180W.c20121120.avg_dnb.tif.gz',sep='')
download.file(
'http://www.worldgrids.org/lib/exe/fetch.php?media=lnmdms3a.tif.gz',
nightFile)
library("R.utils")
gunzip(nightFile,overwrite=TRUE)
nightFull = raster(gsub("\\.gz$", "", nightFile))

border2 = spTransform(torontoBorder, CRS(projection(nightFull)))
toMap2 = openmap(border2)

torontoNight = crop(nightFull, raster::extend(extent(border2),0.1))
torontoNight = projectRaster(torontoNight, crs=CRS(proj4string(murder)))
cscaleNight = colourScale(torontoNight,breaks=9,style='equal',dec=0)

map.new(torontoBorder)
plot(toMap,add=TRUE)
plot(torontoNight,add=TRUE, 
	col=cscaleNight$col,
	breaks=cscaleNight$breaks, 
	legend=FALSE,
	alpha=0.5)
plot(torontoBorder,add=TRUE)
legendBreaks('bottomright', cscaleNight)

save(torontoIncome, torontoNight, 
	torontoPdens,
	file="../pkg/geostatsp/data/torontoPop.RData",
	compress='xz')

} # end don't run

}




\keyword{datasets}
