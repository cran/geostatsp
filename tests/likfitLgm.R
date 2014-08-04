library('geostatsp')
n=100
mydat = SpatialPointsDataFrame(cbind(runif(n), runif(n)), 
		data=data.frame(cov1 = rnorm(n), cov2 = rpois(n, 0.5))
)

# get rid of points too close together
thedist =spDists(mydat)
thedist[lower.tri(thedist, diag=TRUE)] = NA
thedist = thedist < 0.01
thedist = apply(thedist, 1, any, na.rm=T)
mydat = mydat[!thedist,]

	

trueParamAniso = param=c(variance=2^2, range=0.2, shape=2,
		nugget=1^2,anisoRatio=4,anisoAngleDegrees=10, nugget=0)



mydat$U = geostatsp::RFsimulate(trueParamAniso,mydat)$variable1
mydat$Y = -3 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
		mydat$U + rnorm(length(mydat), 0, sd=sqrt(trueParamAniso["nugget"]))

mydat$Ybc = (mydat$Y*0.5+1)^2

 

myres = likfitLgm(Ybc ~ cov1 + cov2, mydat, 
		param=c(range=0.1,nugget=0,shape=2, 
				anisoAngleDegrees=20, anisoRatio=2,
				boxcox=0.4), 
		paramToEstimate = c("range","nugget",
				"anisoRatio","anisoAngleDegrees",
				"boxcox") 
)

myres$summary

pdf("ligfitLgm.pdf")
par(mfrow=c(1,2))

myraster = raster(nrows=30,ncols=30,xmn=0,xmx=1,ymn=0,ymx=1)
covEst = matern(myraster, c(0.5, 0.5), par=myres$param)
covTrue = matern(myraster, c(0.5, 0.5), par=trueParamAniso)

plot(covEst, main="estimate")
plot(covTrue, main="true")

dev.off()

library("geostatsp")
data("swissRain")


sr2 = swissRain
sr2$elev = raster::extract(swissAltitude, sr2)
swissFitAgain = likfitLgm(data=sr2, 
		formula=rain~ elev,
		param=c(range=1000,shape=1,nugget=0,boxcox=0.5),
		paramToEstimate = c("range","nugget")
)
swissFitAgain$par		


# test parallel

n=500
mydat = SpatialPointsDataFrame(cbind(runif(n), runif(n)), 
		data=data.frame(cov1 = rnorm(n), cov2 = rpois(n, 0.5))
)

# simulate a random field
trueParam = c(variance=2^2, range=0.15, shape=2, nugget=0.5^2)

mydat$U = geostatsp::RFsimulate(trueParam,mydat)$variable1

# add fixed effects
mydat$Y = -3 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
		mydat$U + rnorm(length(mydat), 0, sd=sqrt(trueParam["nugget"]))


unix.time( likfitLgm(Y ~ cov1 + cov2, mydat, 
				param=c(range=0.1,nugget=0.1,shape=2), 
				paramToEstimate = c("range","nugget")
		)
)

if(FALSE) { 
options(mc.cores = 1)

unix.time( likfitLgmG( Y ~ cov1 + cov2, mydat,
				param=c(range=0.1,nugget=0.1,shape=2), 
				paramToEstimate = c("range","nugget")
		)
)

options(mc.cores = 2)

unix.time( likfitLgmG( Y ~ cov1 + cov2, mydat,
				param=c(range=0.1,nugget=0.1,shape=2), 
				paramToEstimate = c("range","nugget")
		)
)
}