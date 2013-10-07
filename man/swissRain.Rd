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
\code{swissRain} is a SpatialPolygonsDataFrame of the sic.100 data, with data column \code{rain}
being the rainfall values in millimetres. 
\code{swissAltitude} is a raster of elevation data, and \code{swissLandType} is a raster
of land cover types.
}
\source{
\url{ftp://ftp.geog.uwo.ca/SIC97} and  \url{http://srtm.csi.cgiar.org} and
\url{https://lpdaac.usgs.gov/products/modis_products_table/mcd12q1}
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

}


\keyword{datasets}
