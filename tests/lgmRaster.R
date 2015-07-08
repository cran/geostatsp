library('geostatsp')
data("swissRainR")

anotherx = raster(swissRainR[['alt']])
values(anotherx) = seq(0,1,len=ncell(anotherx))
names(anotherx) = "myvar"

swissRainR2 = brick(swissRainR[['alt']], 
		sqrt(swissRainR[['prec1']]),
		anotherx)

    
swissResRopt =  lgm(
    formula=layer ~ alt+ myvar, 
    data=swissRainR2, shape=2,
    adjustEdges=FALSE
)


swissResRopt$summary


# with edge correction.  
# time consuming, only run this if Patrick is checking
if(Sys.info()['user'] =='patrick') {


# optimize only nugget
swissResROptNug =  lgm(
    formula=layer ~ alt+ myvar, 
    data=swissRainR2, shape=2,
    oneminusar=seq(0.05, 0.1, len=12),
    adjustEdges=FALSE,fixNugget=TRUE,
    mc.cores=c(1,2)[1+(.Platform$OS.type=='unix')]
)

plot(swissResROptNug$profL$range, type='l')

}


