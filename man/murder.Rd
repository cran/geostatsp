\name{murder}
\alias{murder}
\alias{torontoPdens}
\alias{torontoIncome}
\alias{torontoNight}
\docType{data}
\title{
Murder locations
 }
\description{
  Locations of murders in Toronto 2005-2010
}
\usage{data(murder)}
\format{
  \code{murder} is a  SpatialPoints object of murder locations.  \code{torontoPdens}, 
  \code{torontoIncome}, and \code{torontoNight} are rasters containing
  population density (per hectare), median household income, and ambient light 
  respectively.

}

\source{
www3.thestar.com/static/googlemaps/homicidemap.html 
}

\examples{
data(murder)
plot(torontoIncome)
points(murder, col="#0000FF40", cex=0.5)

\dontrun{
library(XML)
murder = xmlTreeParse("http://www3.thestar.com/static/googlemaps/homicides.xml" )
murder = as.data.frame(t(xmlSApply(murder, function(qq) as.numeric(xmlAttrs(qq[["point"]])))))
names(murder) = c("long","lat")
murder= na.omit(murder)
murder = SpatialPoints(murder, proj4string = CRS("+init=epsg:4326"))
}


}




\keyword{datasets}
