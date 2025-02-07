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
\code{swissRain} is a \code{SpatVector} 
100 daily rainfall
measurements made in Switzerland on the 8th of May 1986.
\code{swissAltitude} is a raster of elevation data, and \code{swissLandType} is a raster
of land cover types.
}
\source{
\url{https://web.archive.org/web/20241008015622/https://wiki.52north.org/AI_GEOSTATS/AI_GEOSTATSData} and  \url{https://srtm.csi.cgiar.org} and
\url{https://web.archive.org/web/20241129120557/https://lpdaac.usgs.gov/data/}
}

\examples{
data("swissRain")
swissRain = unwrap(swissRain)
swissAltitude = unwrap(swissAltitude)
swissBorder = unwrap(swissBorder)
swissLandType = unwrap(swissLandType)
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


legend("left",fill=thecol[names(commonValues)],
		legend=substr(levels(swissLandType)[[1]][
						match(as.integer(names(commonValues)),
								levels(swissLandType)[[1]]$ID),
						"Category"], 1,12),
						bg= 'white'
				)

}


\keyword{datasets}
