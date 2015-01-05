library('geostatsp')
mymodel = c(mean=-0.5, variance=1, 
				range=2, shape=2)

myraster = raster(nrows=15,ncols=20,xmn=0,xmx=10,ymn=0,ymx=7.5)

# some covariates, deliberately with a different resolution than myraster
covA = covB = myoffset = raster(extent(myraster), 10, 10)
values(covA) = as.vector(matrix(1:10, 10, 10))
values(covB) = as.vector(matrix(1:10, 10, 10, byrow=TRUE))
values(myoffset) = round(seq(-1, 1, len=ncell(myoffset)))

myCovariate = list(a=covA, b=covB, offsetFooBar = myoffset)

set.seed(0)
myLgcp=simLgcp(mymodel, myCovariate, betas=c(a=-0.1, b=0.25), 
	offset='offsetFooBar',
	rasterTemplate=myraster)

if(require("INLA", quietly=TRUE)) {
res = lgcp(data=myLgcp$events, 
		formula = ~ a + b + offset(offsetFooBar),
		grid=40, 
		covariates=myCovariate,
		buffer=1,
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41)),
		control.mode=list(theta=c(0.18, 0.52),restart=TRUE)
)

res$parameters$summary[,c(1,3,5)]

lgcpRoc =  spatialRoc(res, 
	rr=c(1,1.2, 1.5,2), 
	truth=myLgcp, 
	random=FALSE)
	
dimnames(lgcpRoc)[-1]	
	
plot(lgcpRoc[,1,'onemspec'] , 
	lgcpRoc[,1,'sens'], 
	type='l', 
	xlim=c(0,1), ylim=c(0,1),
	ylab='sensitivity', xlab='1-specificity')
	


}

 





	