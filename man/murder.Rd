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
  Locations of murders in Toronto 1990-2014
}
\usage{data("murder")}
\format{
  \code{murder} is a  \code{SpatVector} object of murder locations.  \code{torontoPdens}, 
  \code{torontoIncome}, and \code{torontoNight} are rasters containing
  population density (per hectare), median household income, and ambient light 
  respectively.  \code{torontoBorder} is a \code{SpatVector}  of the boundary of
  the city of Toronto.

}

\source{
Murder data:\url{https://geo.scholarsportal.info/#r/details/_uri@=1882586248}, 

Lights: \url{https://www.ngdc.noaa.gov/eog/viirs/download_ut_mos.html}

Boundary files: \url{https://www150.statcan.gc.ca/n1/en/catalogue/92-160-X}

Income: \url{https://www150.statcan.gc.ca/n1/en/catalogue/97-551-X2006007}
}

\examples{
data("murder")
murder= unwrap(murder)
torontoBorder = unwrap(torontoBorder)

plot(torontoBorder)
points(murder, col="#0000FF40", cex=0.5)

data("torontoPop")
torontoNight = unwrap(torontoNight)
torontoIncome = unwrap(torontoIncome)
torontoPdens = unwrap(torontoPdens)

# light
plot(torontoNight, main="Toronto ambient light")
plot(torontoBorder, add=TRUE)
points(murder, col="#0000FF40", cex=0.5)


# income
plot(torontoIncome, main="Toronto Income")
points(murder, col="#0000FF40", cex=0.5)
plot(torontoBorder, add=TRUE)

# population density
plot(torontoPdens, main="Toronto pop dens")
points(murder, col="#0000FF40", cex=0.5)
plot(torontoBorder, add=TRUE)

}




\keyword{datasets}
