# number of cells... smaller is faster but less interesting
Ncell = 25

# as in example
require('geostatsp')
 
data('swissRain')
swissRain$lograin = log(swissRain$rain)

if(require('INLA', quietly=TRUE)) {
swissFit =  glgm("lograin", swissRain, Ncell, 
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", param=c(.1, .1))))
)

swissFit$parameters$summary
pdf("swissGlgmExc.pdf")
swissExc = excProb(swissFit$inla$marginals.random$space, 0, template=swissFit$raster)
plot(swissExc, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
		col=c('green','yellow','orange','red'))	
plot(swissBorder, add=TRUE)		
dev.off()
pdf("swissGlgmExc2.pdf")
swissExcP = excProb(swissFit$inla$marginals.predict, 3, template=swissFit$raster)
plot(swissExcP, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
		col=c('green','yellow','orange','red'))	
plot(swissBorder, add=TRUE)		
dev.off()
# intercept only
swissFit =  glgm(lograin~1,swissRain, Ncell, 
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", param=c(.1, .1))))
)

swissFit$parameters$summary


	swissExc = excProb(swissFit$inla$marginals.random$space, 0, template=swissFit$raster)
	plot(swissExc, breaks = c(0, 0.2, 0.8, 0.95, 1.00001), 
		col=c('green','yellow','orange','red'))	
	plot(swissBorder, add=TRUE)		
 


# now with formula
swissFit =  glgm(lograin~ CHE_alt,
		swissRain, 
		Ncell, 
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)
swissFit$parameters$summary

# formula, named list elements
swissFit =  glgm(lograin~ elev,
		swissRain, Ncell, 
		covariates=list(elev=swissAltitude), 
		family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)
swissFit$parameters$summary

# categorical covariates
swissFit =  glgm(lograin ~ elev + factor(land),
		swissRain, Ncell, 
covariates=list(elev=swissAltitude,land=swissLandType), 
family="gaussian", buffer=20000,
priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
control.family=list(hyper=list(prec=list(prior="loggamma", 
						param=c(.1, .1))))
)
swissFit$parameters$summary
table(swissFit$inla$.args$data$land)
plot(swissFit$raster[['predict.mean']])

# put some missing values in covaritates
# also don't put factor() in formula
temp = values(swissAltitude)
temp[seq(10000,12000)] = NA
values(swissAltitude) = temp
swissFit =  glgm(rain ~ elev + land,swissRain,  Ncell, 
		covariates=list(elev=swissAltitude,land=swissLandType), 
		family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.2, 2), range=c(50000,500000)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE),
		control.family=list(hyper=list(prec=list(prior="loggamma", 
								param=c(.1, .1))))
)
swissFit$parameters$summary



data('loaloa')
rcl = rbind(
		# wedlands and mixed forests to forest
		c(5,2),c(11,2),
# savannas to woody savannas
		c(9,8),
		# croplands and urban changed to crop/natural mosaid
		c(12,14),c(13,14))
ltLoaR = reclassify(ltLoa, rcl)
levels(ltLoaR) = levels(ltLoa)

 
elevationLoa = elevationLoa - 750
elevLow = reclassify(elevationLoa, c(0, Inf, 0))
elevHigh = reclassify(elevationLoa, c(-Inf, 0, 0))

 covList = list(elLow = elevLow, elHigh = elevHigh, 
		land = ltLoaR, evi=eviLoa)

 loaFit = glgm(
		 y ~ land + evi + elHigh + elLow, #+ f(villageID,model="iid"),
		 loaloa,
		  Ncell, 
		  covariates=covList, 
		family="binomial", Ntrials = loaloa$N,
		shape=2, buffer=25000,
		priorCI = list(sd=c(0.2, 4), range=c(20000,500000)))

loaFit$par$summary

png("loaFitted.png")
plot(loaFit$raster[["predict.invlogit"]])
dev.off()

# prior for observation standard deviation
swissFit =  glgm( formula="lograin",data=swissRain, grid=Ncell,
		covariates=swissAltitude, family="gaussian", buffer=20000,
		priorCI=list(sd=c(0.1, 2), range=c(50000,500000), 
				sdNugget=c(0.1, 2)), 
		control.mode=list(theta=c(1.9,0.15,2.6),restart=TRUE)
)



# a model with little data, posterior should be same as prior

data2 = SpatialPointsDataFrame(cbind(c(1,0), c(0,1)),
		data=data.frame(y=c(0,0), offset=c(-100,-200)))

res = glgm(data=data2, grid=20, formula=y~1, 
covariates=NULL,
		priorCI = list(sd=c(1,2), range=c(0.3, 2)),
		family="poisson",buffer=0.5)

priorPrec = res$par$sd$params
priorRange = res$par$range$params
pdf("nodata.pdf")

par(mfrow=c(2,1))
# sd
plot(res$parameters$sd$prior,type='l', col='blue')
lines(res$parameters$sd$post,col='red')


# range
plot(res$parameters$range$prior,type='l', col='blue')
lines(res$parameters$range$post,col='red')
legend("topright", col=c("blue","red"),lty=1,legend=c("prior","post'r"))
dev.off()

 

}
