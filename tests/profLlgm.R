library('geostatsp')
data('swissRain')

Ncores = c(1,2)[1+(.Platform$OS.type=='unix')]

swissFit = lgm(data=swissRain, formula=rain~ SRTM_1km,
		grid=150, covariates=swissAltitude,
		shape=1,  fixShape=TRUE, 
		boxcox=0.5, fixBoxcox=TRUE, 
		aniso=TRUE,reml=TRUE,
		param=c(anisoAngleDegrees=37,anisoRatio=10))


x=profLlgm(swissFit, mc.cores=Ncores,
		anisoAngleDegrees=seq(30, 43 , len=12)
)


swissInf = informationLgm(swissFit)



pdf("profLswissAngle.pdf")

plot(x[[1]],x[[2]], xlab=names(x)[1],
#		yaxt='n',
		ylab='log L',
		ylim=c(min(x[[2]]),x$maxLogL),
		type='n')
abline(h=x$breaks[-1],
		col=x$col,
		lwd=1.5)
axis(2,at=x$breaks,labels=x$prob,line=-1.2,tick=F,
		las=1,padj=1.2,hadj=0,col.axis='red')

abline(v=x$ciLong$par,
		lty=2,
		col=x$col[as.character(x$ciLong$prob)])


axis(1,at=x$ciLong$par,
		labels=x$ciLong$quantile,
		padj= -6,hadj=0.5, 
		tcl=0.5,cex.axis=0.8,
		col=NA,col.ticks='red',col.axis='red')

ciCols = grep("^ci", colnames(swissInf$summary),
		value=TRUE)
axis(1,at=swissInf$summary[names(x)[[1]],ciCols],
		labels=gsub("^ci","",ciCols),
		padj= 2,hadj=0.5, 
		tcl=-2,cex.axis=0.7,
		col=NA,col.ticks='blue',col.axis='blue')

lines(x[[1]],x[[2]])


dev.off()

x2d=profLlgm(swissFit, mc.cores=Ncores,
		anisoAngleDegrees=seq(30, 43 , len=6),
		anisoRatio = exp(seq(log(3.5),log(18),len=8))
)
pdf("profLswiss2d.pdf")
image(x2d[[1]],x2d[[2]],x2d[[3]],
		breaks=x2d$breaks,
		col=x2d$col,log='y',
		xlab=names(x2d)[1],
		ylab=names(x2d)[2])

thesevars = c("anisoAngleDegrees","log(anisoRatio)")
thisV = swissInf$information[
		thesevars,thesevars]
thisMean= c(x2d$MLE["anisoAngleDegrees"],
		log(x2d$MLE['anisoRatio']))

haveEllipse = 'ellipse' %in% installed.packages()[,'Package']
if(haveEllipse) {
library('ellipse')

for(D in x2d$prob[x2d$prob>0&x2d$prob<1]) {
	thisE = ellipse(thisV, centre=thisMean,
			level=D)
	thisE = cbind(thisE,
			anisoRatio = exp(thisE[,"log(anisoRatio)"]))
	lines(thisE[,"anisoAngleDegrees"],
			thisE[,"anisoRatio"],lwd=4)
	lines(thisE[,"anisoAngleDegrees"],
			thisE[,"anisoRatio"], col=x2d$col[as.character(D)],
			lwd=3)
}
}

points(x2d$MLE[1],x2d$MLE[2],pch=15) 


if('mapmisc' %in% installed.packages()[,'Package']) {

library('mapmisc')
legendBreaks("topleft",x2d$prob,
		col=x2d$col)
}



dev.off()