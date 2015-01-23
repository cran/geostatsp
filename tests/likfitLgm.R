library('geostatsp')
n=100

set.seed(0)
mydat = SpatialPointsDataFrame(cbind(seq(0,1,len=n), runif(n)), 
		data=data.frame(cov1 = rnorm(n), cov2 = rpois(n, 0.5))
)

# get rid of points too close together
thedist =spDists(mydat)
thedist[lower.tri(thedist, diag=TRUE)] = NA
thedist = thedist < 0.01
thedist = apply(thedist, 1, any, na.rm=T)
mydat = mydat[!thedist,]

	

trueParamAniso = param=c(variance=2^2, range=0.2, shape=2,
		nugget=1^2,anisoRatio=4,anisoAngleDegrees=10)



mydat$U = RFsimulate(trueParamAniso,mydat)$sim

mydat$Y = -3 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
		mydat$U + rnorm(length(mydat), 0, sd=sqrt(trueParamAniso["nugget"]))

mydat$Ybc = (mydat$Y*0.5+1)^2

 
print(range(mydat$Ybc))

date()
myres = likfitLgm(Ybc ~ cov1 + cov2, mydat, 
		param=c(range=0.1,nugget=0,shape=2, 
				anisoAngleDegrees=20, anisoRatio=2,
				boxcox=0.4), 
		paramToEstimate = c("range","nugget",
				"anisoRatio","anisoAngleDegrees",
				"boxcox","shape") 
)
date()

myres$summary[,grep("^ci", colnames(myres$summary),invert=TRUE)]

loglikLgm(formula=Ybc ~ cov1 + cov2, 
    data=mydat, 
    param=myres$param
)


pdf("ligfitLgm.pdf")
par(mfrow=c(1,2))

myraster = raster(nrows=30,ncols=30,xmn=0,xmx=1,ymn=0,ymx=1)
covEst = matern(myraster, y=c(0.5, 0.5), par=myres$param)
covTrue = matern(myraster, y=c(0.5, 0.5), par=trueParamAniso)

plot(covEst, main="estimate")
plot(covTrue, main="true")

dev.off()

library("geostatsp")
data("swissRain")


if(interactive()  | Sys.info()['user'] =='patrick') {
sr2 = swissRain
sr2$elev = raster::extract(swissAltitude, sr2)
swissFitAgain = likfitLgm(data=sr2, 
		formula=rain~ elev,
		param=c(range=1000,shape=1,nugget=0.1,boxcox=0.5),
		paramToEstimate = c("range","nugget")
)
swissFitAgain$par		
}