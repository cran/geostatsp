library('geostatsp')
data("swissRainR")

anotherx = raster(swissRainR[['alt']])
values(anotherx) = seq(0,1,len=ncell(anotherx))
names(anotherx) = "myvar"

swissRainR2 = brick(swissRainR[['alt']], 
		sqrt(swissRainR[['prec1']]),
		anotherx)

swissResR =  lgm( formula=layer ~ alt+ myvar, 
		data=swissRainR2, shape=2,
		oneminusar=seq(0.05, 0.1, len=6),
		nugget =  seq(0.0,0.01,len=20),
		adjustEdges=FALSE,
		mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
)

swissResR$summary[c('oneminusar','range','propNugget',
				grep("\\.betaHat$", rownames(swissResR$summary), value=TRUE)),]

# profile likelihood plot

# plot map of predicted values

swissResR =  lgm( formula=layer ~ alt+ myvar,  
		data=swissRainR2, shape=2,
		oneminusar=seq(0.05, 0.1, len=3),
		nugget =  seq(0,0.01,len=5),
		adjustEdges=TRUE,
		mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
)

swissResR$summary[c('oneminusar','range','propNugget',
				grep("\\.betaHat$", rownames(swissResR$summary), value=TRUE)),]

# range in km
swissResR$summary[ 'range' ,] * sqrt(mean(values(area(swissRainR))))/mean(res(swissRainR))

# profile likelihood plot

# plot map of predicted values


# a simulation
if(FALSE) {
	
	myModel = c(intercept=0,variance=2^2,nugget=3^2, range=2.5,shape=2, 
			cov1=0.2, cov2=-0.5)
	data('nn32')
	theNN=nn64
	covariates = brick(
			xmn=0,ymn=0,xmx=10,ymx=10,
			ncols=attributes(theNN)$Nx,
			nrows=attributes(theNN)$Ny,
			nl=2)
	values(covariates)[,1] = rep(seq(0,1,len=nrow(covariates)), ncol(covariates))
	values(covariates)[,2] = rep(seq(0,1,len=nrow(covariates)), 
			rep(nrow(covariates), ncol(covariates)))
	names(covariates) = c("cov1","cov2")
	myU = RFsimulate(model=myModel,x=raster(covariates), n=5)
	myMean = myModel["intercept"] 
	for(D in names(covariates))
		myMean = myMean + covariates[[D]]*myModel[D]
	myY = myU + myMean
	values(myY) = values(myY)+ rnorm(ncell(myY), 0, sqrt(myModel['nugget']))

	
	covMat = as.matrix(cbind(intercept=1,as.data.frame(covariates)))

	res = loglikGmrf(
			oneminusar=seq(0.01, 0.1, len=4),
			Yvec=as.data.frame(myY), 
			Xmat=covMat, 
			NN=theNN, 
			propNugget=seq(0,2,len=5),
			shape=1, 
			adjustEdges=FALSE,
			boxcox=1, reml=FALSE,
			mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
	)
		
	par=summaryGmrfFit(res)[,'mle','ml',]
	
	myCond = conditionalGmrf(
			param=par,
			Yvec=as.data.frame(myY),
			Xmat=covMat, 
			NN=theNN,
			template=raster(covariates)
			)

		plot(myCond[[grep("2_predict$", names(myCond))[1]]])	
			
			# now with only one simulation
			resOne = loglikGmrf(
					oneminusar=seq(0.01, 0.1, len=4),
					Yvec= as.data.frame(myY)[,1], 
					Xmat=covMat, 
					NN=theNN, 
					propNugget=seq(0,2,len=5),
					shape=1, 
					adjustEdges=FALSE,
					boxcox=1, reml=FALSE,
					mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
			)
			
			parOne=summaryGmrfFit(resOne)[,'mle','ml']
			
			myCond = conditionalGmrf(
					param=parOne,
					Yvec=as.data.frame(myY)[,1],
					Xmat=covMat, 
					NN=theNN,
					template=raster(covariates),
					mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
			)
			plot(myCond[['predict']])	
			
}
