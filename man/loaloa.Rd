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
mapmisc::legendBreaks("bottomleft",ltLoa, bty='n')
points(loaloa)

}


\keyword{datasets}
