\name{spatialRoc}
\alias{spatialRoc}

\title{
Sensitivity and specificity
}
\description{
Calculate ROC curves using model fits to simulated spatial data
}
\usage{
spatialRoc(fit, rr = c(1, 1.2, 1.5, 2), truth, border=NULL, random = FALSE)
}
\arguments{
  \item{fit}{
	A fitted model from the \code{\link{lgcp}} function
}
  \item{rr}{
	Vector of relative risks exceedance probabilities will be calculated for.  Values
	are on the natural scale, with \code{spatialRoc} taking logs when appropriate.
}
  \item{truth}{
True value of the spatial survace, or result from \code{\link{simLgcp}} function.
Assumed to be on the log scale if \code{random=TRUE} and on the natural scale otherwise.
}
  \item{border}{
optional, \code{SpatialPolygonsDataFrame} specifying region that calculations will be restricted to.
}
  \item{random}{
compute ROC's for relative intensity (\code{FALSE}) or random effect (\code{TRUE})
}
}
\details{
Fitted models are used to calculate exceedance probabilities, and
a location is judged to be above an \code{rr} threshold if this 
exceedance probability is above a specified probability threshold. 
Each raster cell of the true surface is categorized as being either true positive, false
positive, true negative, and false negative and sensitivity and specificity computed.
ROC curves are produced by varying the probability threshold.  
}
\value{
An array, with dimension 1 being probability threshold, dimension 2 
being the relative risk threshold, dimension 3 being sensitivity and specificity.
If \code{fit} is a list of model fits, dimension 4 corresponds to elements of \code{fit}.
}
\author{
Patrick Brown 
}


\seealso{
\code{\link{lgcp}}, \code{\link{simLgcp}}, \code{\link{excProb}}
}
 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
 