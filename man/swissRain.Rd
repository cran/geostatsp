\name{swissRain}
\alias{swissRain}
\alias{swissAltitude}
\alias{swissBorder}

\docType{data}
\title{
Swiss rainfall data
}
\description{
Data from the SIC-97 project: Spatial Interpolation Comparison.
}
\usage{data(swissRain)}
\format{
\code{swissRain} is a SpatialPolygonsDataFrame of the sic.100 data, with data column \code{rain}
being the rainfall values in millimetres. 
\code{swissAltitude} is a raster of elevation data
}
\source{
\url{ftp://ftp.geog.uwo.ca/SIC97} and  \url{http://srtm.csi.cgiar.org}
}

\examples{
data(swissRain)
plot(swissAltitude)
points(swissRain)
plot(swissBorder, add=TRUE)



}


\keyword{datasets}
