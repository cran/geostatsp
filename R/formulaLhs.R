
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


