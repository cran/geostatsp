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
\dontrun{
library('raster')
border = raster::getData('GADM',country='CHE',level=0)

myextent = extent(5 ,12,45,49)

theres=10
rain = raster::getData('worldclim', var='prec',
	res=theres,mask=FALSE)#,lon=20,lat=47)
rain = raster::crop(rain,myextent)
plot(rain[['prec3']])
plot(border,add=TRUE)

alt = raster::getData('worldclim', var='alt',res=theres)
alt = raster::crop(alt,extent(rain))
plot(alt)
plot(border,add=TRUE)


bio = raster::getData('worldclim', var='bio',res=theres)
bio = raster::crop(bio,extent(rain))
plot(alt)
plot(border,add=TRUE)

swissRainR = rain
#	rain[[c('prec2','prec7')]]
swissRainR = addLayer(swissRainR,alt)	

library('geostatsp')
swissNN = NNmat(swissRainR)
	
save(swissRainR, swissNN,file=
"/home/patrick/workspace/diseasemapping/pkg/gmrf/data/swissRainR.RData",
compress='xz'	)	

}
data('swissRainR')
plot(swissRainR[['prec7']])

plot(swissRainR[['alt']])

swissNN[1:4,1:5]

}
\keyword{datasets}
