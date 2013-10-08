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
\url{http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.005/} for the EVI data,  \url{http://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.051/} for land type and  
\url{http://srtm.csi.cgiar.org} for the elevation data.
}

\examples{


data("loaloa")
plot(loaloa, main="loaloa villages")

# elevation
data("elevationLoa")
plot(elevationLoa, col=terrain.colors(100), main="elevation")
points(loaloa)

# vegetation index
data("eviLoa")
plot(eviLoa, main="evi")
points(loaloa)

data("tempLoa")
plot(tempLoa, main="temperature")
points(loaloa)

data("ltLoa")
# land type, a categorical variable
commonValues  = sort(table(values(ltLoa)),decreasing=TRUE)[1:7]
commonValues=commonValues[!names(commonValues)==0]

thelevels = levels(ltLoa)[[1]]$ID
thebreaks = c(-0.5, 0.5+thelevels)
thecol = rep(NA, length(thelevels))
names(thecol) = as.character(thelevels)

thecol[names(commonValues)] = rainbow(length(commonValues))
plot(ltLoa, breaks=thebreaks, col=thecol,legend=FALSE, main="land type")
legend("bottomleft",fill=thecol[names(commonValues)],
		legend=levels(ltLoa)[[1]][
						match(as.integer(names(commonValues)),
								levels(ltLoa)[[1]]$ID),
						"Category"],
				bty='n'
				)

points(loaloa)



}


\keyword{datasets}
