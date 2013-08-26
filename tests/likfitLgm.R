library(geostatsp)
n=100
mydat = SpatialPointsDataFrame(cbind(runif(n), runif(n)), 
		data=data.frame(cov1 = rnorm(n), cov2 = rpois(n, 0.5))
)


 trueParamAniso = param=c(variance=2^2, range=0.2, rough=2,
		nugget=0,aniso.ratio=4,aniso.angle.degrees=10, nugget=0)

mydat$U = GaussRF(mydat, par=trueParamAniso)
mydat$Y = -3 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
		mydat$U + rnorm(length(mydat), 0, sd=sqrt(trueParamAniso["nugget"]))

mydat$Ybc = (mydat$Y*0.5+1)^2

 

myres = likfitLgm(mydat, Ybc ~ cov1 + cov2, 
		param=c(range=0.1,nugget=0,rough=2, 
				aniso.angle.degrees=0, aniso.ratio=2,
				boxcox=0.4), 
		paramToEstimate = c("range","nugget",
				"aniso.ratio","aniso.angle.degrees",
				"boxcox") 
)

myres$summary

par(mfrow=c(1,2))

myraster = raster(nrows=30,ncols=30,xmn=0,xmx=1,ymn=0,ymx=1)
covEst = matern(myraster, c(0.5, 0.5), par=myres$param)
covTrue = matern(myraster, c(0.5, 0.5), par=trueParamAniso)

plot(covEst, main="estimate")
plot(covTrue, main="true")

par(mfrow=c(1,1))


