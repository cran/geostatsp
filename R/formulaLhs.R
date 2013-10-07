
formulaLhs = function(formula) {

	observations = gsub("^~,[[:space:]]+", "", toString(formula))
	observations = gsub(",[[:print:]]+$", "", observations)
	
	observations
}

formulaRhs = function(formula, char=FALSE) {

	x = gsub("^.*~", "", toString(base::format(formula)))

	if(!char)
		x= as.formula(paste("~", x))
	x
}


