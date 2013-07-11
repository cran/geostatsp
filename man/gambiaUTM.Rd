\name{gambiaUTM}
\alias{gambiaUTM}
\docType{data}
\title{
Gambia data
}
\description{
  This data-set was used by Diggle, Moyeed, Rowlingson, and Thomson (2002) to demonstrate 
  how the model-based geostatistics framework of Diggle et al. (1998) could be adapted to 
  assess the source(s) of extrabinomial variation in the data and, in particular, whether this 
  variation was spatially structured. The malaria prevalence data set consists of measurements 
  of the presence of malarial parasites in blood samples obtained from children in 65 villages 
  in the Gambia. Other child- and village-level indicators include age, bed net use, whether the
  bed net is treated, whether or not the village belonged to the primary health care structure,
  and a measure of 'greenness' using a vegetation index.
}
\usage{data(gambiaUTM)}
\format{
  A SpatialPointsDataFrame, with column \code{pos} being the binary response for a malaria 
  diagnosis, as well as other child-level indicators such as \code{netuse} and \code{treated} 
  being measures of bed net use and whether the nets were treated. The column \code{green} is 
  a village-level measure of greenness. A UTM coordinate reference system is used, where 
  coordinates are in metres.
}

\source{
\url{http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets}.  For further details 
on the malaria data, see Thomson et al. (1999).
}

\examples{
data(gambiaUTM)

plot(gambiaUTM)

\dontrun{
download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS-0.3.zip", 
  "borders.zip")
unzip("borders.zip")
worldBorders = readOGR(".", "TM_WORLD_BORDERS-0.3")
africa = worldBorders[worldBorders$REGION ==2,]
plot(spTransform(africa, gambiaUTM@proj4string), add=T)

}


}

\references{
  Diggle, P. J., Moyeed, R. A., Rowlingson, R. and Thomson, M. (2002). Childhood Malaria in 
  the Gambia: A case-study in model-based geostatistics. Journal of the Royal Statistical 
  Society. Series C (Applied Statistics), 51(4): 493-506.
  
  Diggle, P. J., Tawn, J. A. and Moyeed, R. A. (1998). Model-based
  geostatistics (with Discussion). Applied Statistics,
  47, 299--350.

  Thomson, M. C., Connor, S. J., D'Alessandro, U., Rowlingson, B., Diggle, P., Creswell, 
  M. and Greenwood, B. (2004). Predicting malaria infection in Gambian children from 
  satellite data and bed net use surveys: the importance of spatial correlation in the 
  interpretation of results. American Journal of Tropical Medicine and Hygiene, 61: 2-8.
}


\keyword{datasets}
