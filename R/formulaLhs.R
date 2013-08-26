
formulaLhs = function(formula) {


	observations = gsub("~.*$", "", base::format(formula))
	observations = gsub("^[[:space:]]+|[[:space:]]+$", "", observations)
	
#	observations = base::as.character(formula)
#	observations = observations[-c(1,length(observations))]
	observations
}

formulaRhs = function(formula) {

	x = gsub("^.*~", "~", base::format(formula))
#	x = base::as.character(formula)
#	x=x[length(x)]
#	x = paste("~", x)
	
	as.formula(x)
}