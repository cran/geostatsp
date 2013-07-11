\name{loaloa}
\alias{loaloa}
\alias{elevationLoa}
\alias{eviLoa}
\alias{ltLoa}
\docType{data}
\title{
Loaloa prevalence data from 197 village surveys
}
\description{
Location and prevalence data from villages, elevation an vegetation index for the study region.
}
\usage{data(loaloa)}
\format{
\code{loaloa} is a SpatialPolygonsDataFrame of the data, with columns \code{N} being the number 
of individuals tested and \code{y} being the number of positives.
\code{elevationLoa} is a raster of elevation data.
\code{eviLoa} is a raster of vegetation index for a specific date. \code{ltLoa} is land type.
\code{ltLoa} is a raster of land types. 1     2     5     6     7     8     9    10    11    12    13    14    15
}

\source{
\url{http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets} for the loaloa data,
\url{http://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.005/} for the EVI data,  \url{http://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.051/} for land type and  
\url{http://srtm.csi.cgiar.org} for the elevation data.
}

\examples{
data(loaloa)
plot(ltLoa, breaks = c(-0.5, levels(ltLoa)[[1]]$ID+0.5),
col=levels(ltLoa)[[1]]$rgb, legend=FALSE)

forLegend = as.numeric(names(sort(table(values(ltLoa)),decreasing=TRUE)))[1:4]
forLegend = match(levels(ltLoa)[[1]]$ID , forLegend, nomatch = 0) > 0

legend("topright", fill=levels(ltLoa)[[1]][forLegend,"rgb"], 
	legend=
	substr(as.character(levels(ltLoa)[[1]][forLegend,"CLASSNAMES"]),1,24))

points(loaloa)



}


\keyword{datasets}
