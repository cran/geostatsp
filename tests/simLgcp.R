
if(requireNamespace("INLA", quietly=TRUE) ) {
  INLA::inla.setOption(num.threads=2)
  # not all versions of INLA support blas.num.threads
  try(INLA::inla.setOption(blas.num.threads=2), silent=TRUE)
} 

library('geostatsp')

# exclude this line to use the RandomFields package
options(useRandomFields = FALSE)

mymodel = c(mean=-1.5, variance=1, 
				range=2, shape=2)

myraster = rast(nrows=15,ncols=15,xmin=0,xmax=10,ymin=0,ymax=10)

# some covariates, deliberately with a different resolution than myraster
covA = covB = myoffset = rast(ext(myraster), 10, 10)
values(covA) = as.vector(matrix(1:10, 10, 10))
values(covB) = as.vector(matrix(1:10, 10, 10, byrow=TRUE))
values(myoffset) = round(seq(-1, 1, len=ncell(myoffset)))

myCovariate = list(a=covA, b=covB, offsetFooBar = myoffset)

set.seed(0)
myLgcp=simLgcp(mymodel, myCovariate, 
    betas=c(a=-0.1, b=0.25), 
	offset='offsetFooBar',
	rasterTemplate=myraster)

if(requireNamespace("INLA", quietly=TRUE)) {
res = lgcp(data=myLgcp$events, 
		formula = ~ a + b + offset(offsetFooBar),
		grid=squareRaster(myoffset, 15), 
		covariates=myCovariate,
		prior=list(sd=0.2, range=0.4))

res$parameters$summary[,c(1,3,5)]

lgcpRoc =  spatialRoc(
  fit=res, 
	rr=c(2,3), 
	truth=myLgcp, 
	random=FALSE)

head(lgcpRoc)

plot(lgcpRoc[,'onemspec'] , 
	lgcpRoc[,'2'], 
	type='o', 
	xlim=c(0,1), ylim=c(0,1),
	ylab='sensitivity', xlab='1-specificity'
)

}





	