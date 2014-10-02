# as in example
require('geostatsp')
myPoints = SpatialPoints(cbind(rbeta(100,2,2), rbeta(100,3,4)))
myPoints@bbox = cbind(c(0,0), c(1,1))

mycov = raster(matrix(rbinom(100, 1, 0.5), 10, 10), 0, 1, 0, 1)
names(mycov)="x1"

if(require("INLA", quietly=TRUE)) {

res = lgcp(data=myPoints, grid=20, covariates=mycov,
		formula=~factor(x1),
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

	
# intercept only

res = lgcp(data=myPoints, grid=20, covariates=mycov,
		formula=~1,
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

	# dodgy formula

res = lgcp(data=myPoints, grid=20, covariates=mycov,
		formula=shouldntBeHere~1,
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

	
# some missing values

temp = values(mycov)
temp[1:4] = NA
values(mycov) = temp

res = lgcp(data=myPoints, grid=20, covariates=mycov,
		formula=~factor(x1),
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

}	