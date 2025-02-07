## ----knitr, include=FALSE, tidy=FALSE-----------------------------------------
require('knitr')
# very basic output that doesnt require extra latex packages

knit_hooks$set(source  =function(x, options) {
	paste0(c('\\begin{verbatim}', x, '\\end{verbatim}', ''),
			collapse = '\n')
})
knit_hooks$set(plot = function (x, options) 
		{
			paste0(knitr::hook_plot_tex(x, options), "\n")
	 })
knit_hooks$set(chunk = function(x, options) {x})

hook_output <- function(x, options) {
	if (knitr:::output_asis(x, options)) return(x)
	paste0('\\begin{verbatim}\n', x, '\\end{verbatim}\n')
}

knit_hooks$set(output = hook_output)
knit_hooks$set(message = hook_output)
knit_hooks$set(warning = hook_output)
knitr::opts_chunk$set(out.width='0.48\\textwidth', 
	fig.align='default', fig.height=6, fig.width=6,
	tidy = FALSE)

## ----packages-----------------------------------------------------------------
library("geostatsp")
data('murder')
data('torontoPop')
murder = unwrap(murder)
torontoBorder = unwrap(torontoBorder)
torontoPdens = unwrap(torontoPdens)
torontoIncome = unwrap(torontoIncome)

## ----inlaOpts, include=FALSE--------------------------------------------------
if(requireNamespace("INLA", quietly=TRUE) ) {
  INLA::inla.setOption(num.threads=2)
  # not all versions of INLA support blas.num.threads
  try(INLA::inla.setOption(blas.num.threads=2), silent=TRUE)
}

## ----Covariates, tidy=FALSE---------------------------------------------------
theCrs = paste0("+proj=omerc +lat_0=43.7117469868935 +lonc=-79.3789787759006",
	" +alpha=-20 +gamma=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
murderT = project(murder, theCrs)
borderT = project(torontoBorder, crs(murderT))
borderC = crop(borderT, ext(-12700, 7000, -7500, 3100))

covList = list(
		pop=torontoPdens,
		inc = log(torontoIncome) )

formulaHere = ~ inc + offset(pop, log=TRUE)

## ----lgcpGamma, tidy=FALSE----------------------------------------------------
resG=lgcp(
	formula = formulaHere, 
	data=murderT, 
	grid=squareRaster(borderC, 30), 
	covariates=covList,
	border=borderC, 
	buffer=2000,
	prior = list(
				sd = c(lower = 0.2, upper = 2), 
				range = c(lower = 2, upper=20)*1000),
	control.inla=list(strategy='gaussian'))

## ----lgcpGammaRes-------------------------------------------------------------
if(!is.null(resG$parameters)) {
	knitr::kable(resG$parameters$summary, digits=3)
}

## ----lgcpPc, tidy=FALSE-------------------------------------------------------
resP=lgcp(formulaHere, data=murderT, 
			grid=squareRaster(borderC, 30),
			covariates=covList,
			border=borderC, buffer=2000,
			prior = list(
				sd = c(u=0.5, alpha=0.05), 
				range = c(u=10*1000, alpha = 0.4)),
			control.inla = list(strategy='gaussian')
	) 

## ----lgcpGammaRespc-----------------------------------------------------------
if(!is.null(resP$parameters)) {
	knitr::kable(resP$parameters$summary, digits=3)
}

## ----lgcpTable, tidy=FALSE----------------------------------------------------
sdSeq = seq(0,4,len=501)
rangeSeq = seq(0,15*1000, len=501)
resT=lgcp(formulaHere, 
			data=murderT, 
			grid=squareRaster(borderC, 30),
			covariates=covList,
			border=borderC, buffer=2000,
			prior = list(
					sd = cbind(sdSeq, dexp(sdSeq, 2)), 
					range = cbind(rangeSeq, dexp(rangeSeq, 1/5000))),
			control.inla = list(strategy='gaussian')
)

## ----lgcpTableRes-------------------------------------------------------------
if(!is.null(resT$parameters)) {
	knitr::kable(resT$parameters$summary, digits=3)
}

## ----priorPost, fig.cap="Priors and posteriors", fig.subcap=c("sd", "range"), echo=FALSE----

if(!is.null(resG$parameters)) {
	par(mar=c(3,3,0,0))
	
	matplot(resG$parameters$sd$posterior[,'x'], 
			resG$parameters$sd$posterior[,c('y', 'prior')],
			type='l', lty=1:2,
			col='red', 
			xlab='sd', ylab='dens')
	matlines(resP$parameters$sd$posterior[,'x'], 
			resP$parameters$sd$posterior[,c('y', 'prior')],
			lty=1:2, col='blue' )
	matlines(resT$parameters$sd$posterior[,'x'], 
			resT$parameters$sd$posterior[,c('y', 'prior')],
			lty=1:2, col='green' )

	legend('topleft', lty=c(1,1,1,1,2), 
			col=c('red','blue', 'green','black','black'), 
			legend=c( 'Gamma','PC','table', 'posterior','prior'))


	matplot(resG$parameters$range$posterior[,'x'], 
			resG$parameters$range$posterior[,c('y', 'prior')],
			type='l', lty=1:2,
			col='red', 
			xlab='range', ylab='dens')
	matlines(resP$parameters$range$posterior[,'x'], 
			resP$parameters$range$posterior[,c('y', 'prior')],
			lty=1:2, col='blue' )
	matlines(resT$parameters$range$posterior[,'x'], 
			resT$parameters$range$posterior[,c('y', 'prior')],
			lty=1:2, col='green' )


} else {
	plot(1:10)
	plot(1:10)
}

## ----maps, fig.cap='Random effects and fitted values', fig.subcap=c('gamma, fitted','pc fitted','gamma random','pc random'), fig.height=3, fig.width=5, echo=FALSE----

if(requireNamespace('mapmisc', quietly=TRUE) & !is.null(resG$raster)) {
	
	thecex=1.2	
	
	colFit = mapmisc::colourScale(resG$raster[['predict.exp']],
			breaks=6, dec=8, style='equal',
			transform='log')
	
	
	mapmisc::map.new(borderC, TRUE)
	plot(resG$raster[['predict.exp']], 
			col=colFit$col,breaks=colFit$breaks,
			legend=FALSE, add=TRUE)
	plot(borderC, add=TRUE)
	points(murderT, col='#00FF0050',cex=0.6)
	mapmisc::legendBreaks('bottomright', colFit, cex=thecex)
	
	
	mapmisc::map.new(borderC, TRUE)
	plot(resP$raster[['predict.exp']], 
			col=colFit$col,breaks=colFit$breaks,
			legend=FALSE, add=TRUE)
	plot(borderC, add=TRUE)
	points(murderT, col='#00FF0050',cex=0.6)
	mapmisc::legendBreaks('bottomright', colFit, cex=thecex)
	
	
	colR = mapmisc::colourScale(resG$raster[['random.mean']],
			breaks=12, dec=0, style='equal')
	
	mapmisc::map.new(borderC, TRUE)
	plot(resG$raster[['random.mean']], 
			col=colR$col,breaks=colR$breaks,
			legend=FALSE, add=TRUE)
	plot(borderC, add=TRUE)
	points(murderT, col='#00FF0050',cex=0.6)
	mapmisc::legendBreaks('bottomright', colR, cex=thecex)
	
	mapmisc::map.new(borderC, TRUE)
	plot(resP$raster[['random.mean']], 
			col=colR$col,breaks=colR$breaks,
			legend=FALSE, add=TRUE)
	plot(borderC, add=TRUE)
	points(murderT, col='#00FF0050',cex=0.6)
	mapmisc::legendBreaks('bottomright', colR, cex=thecex)
	
	
} else {
	plot(1:10)
	plot(1:10)
	plot(1:10)
	plot(1:10)
	
}

