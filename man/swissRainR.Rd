\name{swissRainR}
\alias{swissRainR}
\alias{swissNN}
\docType{data}
\title{
Raster of Swiss rain data
}
\description{
A raster image of Swiss rain and elevation, and a nearest neighbour matrix
corresponding to this raster.
}
\usage{data(swissRainR)}
\format{
\code{swissRainR} is a RasterBrick of Swiss elevation and 
precipitation, and \code{swissNN} is a matrix of nearest neighbours.
}

\source{
See examples
}

\examples{

data('swissRainR')
swissRainR = unwrap(swissRainR)
plot(swissRainR[['prec7']])
plot(swissRainR[['alt']])

swissNN[1:4,1:5]

}
\keyword{datasets}
