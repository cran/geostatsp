
formulaLhs = function(formula) {

	observations = gsub("^~,[[:space:]]+", "", toString(formula))
	observations = gsub(",[[:print:]]+$", "", observations)
	
	observations
}

formulaRhs = function(formula, char=FALSE) {
	formula = base::format(formula)
	# if there is a line break in the formula, 
	# format(formula) will create a vector
	formula = paste(formula, collapse="")
	
	
	x = gsub("^.*~", "", toString(formula))


	if(!char)
		x= as.formula(paste("~", x))
	x
}


setGeneric('resampleMethods', function(
				formula, covariates, data=NULL) 
			standardGeneric("resampleMethods"))



setMethod("resampleMethods", 
		signature("formula"), 
		function(formula, covariates, data=NULL){
# decide which method to use when reprojecting covariates
			# factors must be ngb, numerics are bilinear
			allterms = formulaRhs(formula,char=TRUE)
			allterms = unlist(strsplit(allterms, "\\+"))
			allterms = gsub("[[:space:]]", "", allterms)
			allterms = allterms[allterms != "1"]
			
			
			# get rid of offset
			theOffset = grep("^offset\\(", allterms)
			if(length(theOffset)) allterms = allterms[-theOffset]
# get rid of other random effects
			theOtherRE = grep("^f\\(", allterms)
			if(length(theOtherRE)) allterms = allterms[-theOtherRE]
			alltermsWithF = gsub("\\)$", "", allterms)
			allterms = gsub("^factor\\(", "", alltermsWithF)
			
			theFactors = grep("^factor", alltermsWithF, value=T)
			theFactors = gsub("^factor\\(", "", theFactors)
			
# check to see if any rasters are stored as factors
			notFactors = allterms[!allterms %in% theFactors]
			if(length(notFactors)) {
				
				moreFactors = callGeneric(notFactors, covariates)
				
				if(length(moreFactors))
					theFactors = unique(c(theFactors, moreFactors))
			}
			method = rep("bilinear", length(names(covariates)))
			names(method)=names(covariates)
			method[names(method) %in% theFactors] = "ngb" 
			
			method
		}
)

setMethod("resampleMethods", 
		signature("character"), 
		function(formula, covariates, data=NULL){

			# to do: look for factors in data
			
			theFactors = NULL
			formula = match(formula, names(covariates))
			formula = formula[!is.na(formula)]
			for(D in formula) {
				if(any(slotNames(covariates[[D]])=="data")) {
					if(any(slotNames(covariates[[D]]@data)=="isfactor")) {
						if(covariates[[D]]@data@isfactor)					
							theFactors = c(theFactors, D)
					}
				}
			}
			theFactors
		}
)
			