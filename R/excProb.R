excFunQQ = function(themat, threshold) {

	if(length(themat)) {
	over = themat[,"x"]>threshold
	
	if(all(over)){
		prob = 1
	} else if(any(over)) {
		
		toInt = rbind(c(threshold, approx(themat[,"x"], themat[,"y"], threshold)$y),
				themat[over,]
		)
		
		prob = pracma::trapz(toInt[,"x"], toInt[,"y"])
	} else {
		prob = 0
	} 
	} else { # matrix has length zero
		prob = NA
	}
	
	prob
} 

excProb = function(x, threshold=0, random=FALSE, template=NULL, 
		templateIdCol=NULL, nuggetInPrediction=TRUE) {
	

	
elementsColumnwise = FALSE # default is indexes for x[[index]] refer to 
# cells numbered columnwise
	lowertail=FALSE	
	
if(is.list(x))	{
 # model is from lgm
	if(all(c("predict","param", "summary") %in% names(x))){
		template = raster(x$predict)
		if(!random) {
			# check for boxcox
			lowertail=FALSE	
			elementsColumnwise = FALSE
			if(any(names(x$predict)=="predict.boxcox")) {
				resp = x$predict[["predict.boxcox"]]	
				if(x$param["boxcox"]!=0) {
					threshold = (threshold^x$param["boxcox"] - 1) /
						x$param["boxcox"]
					lowertail = x$param["boxcox"] < 0
				} else {
					threshold = log(threshold)
				}
			} else if(any(names(x$predict)=="predict.log")){ 
				resp = log(x$predict[["predict"]])
				threshold = log(threshold)
			} else {
				resp = x$predict[["predict"]]
			}
		}	else { # end if !random
			resp = x$predict[["random"]]
		}
		
		theSD = values(x$predict[["krigeSd"]])
		if(nuggetInPrediction) {
			theSD = sqrt(theSD^2 + x$param["nugget"])
		} 
		
		
		excProbAll = pnorm( threshold , mean=values(resp), sd= theSD,
				lower.tail=lowertail)

		
	} else { # isn't from lgm but is a list
	
	# check if it's from glgm
	 if(all(c("inla", "raster","parameters") %in% names(x))) {
		 template = raster(x$raster)
		if(!random ) {
			x = x$inla$marginals.lincomb.derived
		} else {
			x = x$inla$marginals.random$space
		}
	
	 }
	 excProbAll = unlist(lapply(x, excFunQQ, threshold=threshold))
	 names(excProbAll) = names(x)
	 } # end not from lgm	
		
	} else { # not a list, must be a matrix or data frame with two columns.
		excProbAll = excFunQQ(x, threshold)
		
	}		


	# make sure probabilities are between zero and 1
excProbAll = pmax(0, pmin(excProbAll, 1))

if(length(grep("^Raster", class(template)))) {
	template = raster(template)
	if(elementsColumnwise) {
		values(template) = matrix(excProbAll, 
							nrow=nrow(template),ncol=ncol(template),
							byrow=FALSE)
	} else {
		values(template) = excProbAll
	}
	excProbAll = template
	names(excProbAll) = paste("exc","random"[random], threshold, sep="")
	
} 

if(length(grep("(SpatialPolygonsDataFrame|SpatialPointsDataFrame)", class(template)))) {
	newcol=paste("exc", "random"[random],threshold, sep="")
	if(is.null(templateIdCol)) {
		template[[newcol]] = excProbAll
	} else {
	template[[newcol]] = 
			excProbAll[template@data[templateIdCol]]
	}
	excProbAll = template
} 


excProbAll
} 