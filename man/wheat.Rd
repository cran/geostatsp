\name{wheat}
\alias{wheat}
 
\docType{data}
\title{
Mercer and Hall wheat yield data
}
\description{
Mercer and Hall wheat yield data, based on version in Cressie (1993), p. 455.
}
\usage{data(wheat)}
\format{
\code{wheat} is a raster where the values refer to wheat yields.  
}
 

\examples{
	data("wheat")
	wheat = unwrap(wheat)
	plot(wheat, main="Mercer and Hall Data")
}





\keyword{datasets}
