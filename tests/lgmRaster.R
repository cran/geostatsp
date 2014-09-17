library('geostatsp')
data("swissRainR")
swissRainR2 = brick(swissRainR[['alt']], 
		sqrt(swissRainR[['prec1']]))

swissResR =  lgm( formula=layer ~ alt, 
		data=swissRainR2, shape=2,
		oneminusar=seq(0.05, 0.1, len=6),
		nugget =  seq(0,0.01,len=20),
		adjustEdges=FALSE,
		mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
)

swissResR$summary[c('oneminusar','range','propNugget','(Intercept).betaHat', 'x.betaHat'),]

swissResR =  lgm( formula=layer ~ alt, 
		data=swissRainR2, shape=2,
		oneminusar=seq(0.05, 0.1, len=3),
		nugget =  seq(0,0.01,len=5),
		adjustEdges=TRUE,
		mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
)

swissResR$summary[c('oneminusar','range','propNugget','(Intercept).betaHat', 'x.betaHat'),]

# range in km
swissResR$summary[ 'range' ,] * sqrt(mean(values(area(swissRainR))))/mean(res(swissRainR))