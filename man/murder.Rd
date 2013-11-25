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
  Locations of murders in Toronto 2005-2010
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
http://www.thestar.com/news/crime/torontohomicidemap.html, 
http://ngdc.noaa.gov/eog/viirs/download_viirs_ntl.html,
http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/bound-limit-2011-eng.cfm,
http://www12.statcan.gc.ca/census-recensement/2011/dp-pd/tbt-tt/Index-eng.cfm
}

\examples{
data("murder")
plot(torontoNight, main="Toronto ambient light")
plot(torontoBorder, add=TRUE)
points(murder, col="#0000FF40", cex=0.5)

murder2 = murder[!duplicated(coordinates(murder)),]


data("torontoPop")
plot(torontoIncome, main="Toronto Income")
points(murder2, col="#0000FF40", cex=0.5)

plot(torontoPdens, main="Toronto pop dens")
points(murder2, col="#0000FF40", cex=0.5)


}




\keyword{datasets}
