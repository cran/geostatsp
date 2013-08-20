
excProb = function(marginals, threshold=0, template=NULL, templateIdCol=NULL) {
	

	
excFunQQ = function(themat) {
	if(length(themat)) {
	over = themat[,"x"]>threshold
	
	toInt = rbind(c(threshold, approx(themat[,"x"], themat[,"y"], threshold)$x),
			themat[over,]
			)
	
	prob = trapz(toInt[,"x"], toInt[,"y"])
} else {
	prob = NA
}
	
	prob
} 

if(is.list(marginals)) {
 excProbAll = unlist(lapply(marginals, excFunQQ))
} else {
	excProbAll = excFunQQ(marginals)
}
# make sure probabilities are between zero and 1
excProbAll = pmax(0, pmin(excProbAll, 1))

if(length(grep("^Raster", class(template)))) {
	if(nlayers(template)>1)
		template = template[[1]]
	names(excProbAll) = names(marginals)
	values(template) = matrix(excProbAll, 
							nrow=nrow(template),ncol=ncol(template),
							byrow=T)
	excProbAll = template
	names(excProbAll) = paste("exc", threshold, sep="")
} 

if(length(grep("(SpatialPolygonsDataFrame|SpatialPointsDataFrame)", class(template)))) {
	newcol=paste("exc", threshold, sep="")
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