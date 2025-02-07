
havePackages = c(
    'INLA' = requireNamespace('INLA', quietly=TRUE)
)

if(requireNamespace("INLA", quietly=TRUE) ) {
  INLA::inla.setOption(num.threads=2)
  # not all versions of INLA support blas.num.threads
  try(INLA::inla.setOption(blas.num.threads=2), silent=TRUE)
}

print(havePackages)

# as in example
require('geostatsp')
myPoints = vect(cbind(rbeta(100,2,2), rbeta(100,3,4)))

mycov = rast(matrix(rbinom(100, 1, 0.5), 10, 10), extent=ext(0, 1, 0, 1))
names(mycov)="x1"

res = lgcp(
    formula=~factor(x1),
    data=myPoints, 
    grid=20, 
    covariates=mycov,
		prior=list(sd=0.1, range=c(0.4)),  
)

if(length(res$parameters)) {
knitr::kable(res$parameters$summary, digits=3)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)
}

# intercept only

res = lgcp(
    data=myPoints, 
    grid=20, 
    covariates=mycov,
		formula=~1,
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

if(length(res$parameters)) {

knitr::kable(res$parameters$summary, digits=3)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)
}


# some missing values
mycov2 = deepcopy(mycov)
temp = values(mycov2)
temp[1:4] = NA
values(mycov2) = temp

if(all(havePackages)) {

res = lgcp(data=myPoints, grid=20, covariates=mycov2,
		formula=~factor(x1),
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

if(length(res$parameters)) {

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)
}
}


data('murder')
murder = unwrap(murder)
data('torontoPop')
torontoPdens = unwrap(torontoPdens)
torontoIncome = unwrap(torontoIncome)
torontoBorder = unwrap(torontoBorder)

highDens = torontoPdens > 100


myCov = list(
    pop=torontoPdens,
    highDens = highDens,
    inc = log(torontoIncome)
)

formula = ~ inc*highDens + offset(pop, log=TRUE)


resL=lgcp(formula, data=murder, 
    grid=squareRaster(murder, 30),
    covariates=myCov, verbose=TRUE,
    border=torontoBorder)

resO = lgcp( ~ inc + pop, 
    data=murder, 
    grid=squareRaster(murder, 30),
    covariates=list(inc=myCov$inc, pop=log(myCov$pop)),
    border=torontoBorder)

if(length(resL$parameters)) {

	rbind(resL$param$summary, resO$param$summary)
}



# check spdfToBrock

if(requireNamespace('diseasemapping', quietly=TRUE)){
	require('diseasemapping')
	
	data('kentucky')
	kentucky = unwrap(kentucky)
	
	popList = list(
			'2002' = kentucky[,c('M.50', 'M.55', 'M.60')],
			'2005' = kentucky[,c('F.50', 'F.55', 'F.60')]
			)
	for(D in names(popList))
		names(popList[[D]]) = paste('expected_', seq(as.numeric(D)-1, len=3), sep='')
	
	popBrick = spdfToBrick(
			x=popList,
			template=squareRaster(kentucky, 10),
    	logSumExpected=FALSE
			)
	sum(popList[['2002']]$expected_2001, na.rm=TRUE)		
	sum(values(popBrick[['expected_2001']]), na.rm=TRUE)*prod(res(popBrick))
	
	popBrick2 = spdfToBrick(
			x=popList,
			template=squareRaster(kentucky, 10),
    	logSumExpected=TRUE
	)
	
	sum(unlist(lapply(popList, function(qq) apply(values(qq), 2, sum, na.rm=TRUE))))
	sum(exp(values(popBrick2)), na.rm=TRUE)*prod(res(popBrick2))
			
}
	


