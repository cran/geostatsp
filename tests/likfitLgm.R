library('geostatsp')

# exclude this line to use the RandomFields package
options(useRandomFields = FALSE)

if(interactive()  | Sys.info()['user'] =='patrick') {
  n=300
} else {
  n=100
}

set.seed(0)
mydat = vect(cbind(seq(0,1,len=n), runif(n)), 
		atts=data.frame(cov1 = rnorm(n), cov2 = rpois(n, 0.5))
)


	

trueParamAniso = c(
    variance=2^2, range=0.2, shape=2,
		nugget=1^2,anisoRatio=4,anisoAngleDegrees=10)

trueParamIso = c(
    variance=2^2, range=0.2, shape=2,
    nugget=1^2)


mydat$U = geostatsp::RFsimulate(trueParamAniso,mydat)$sim
mydat$Uiso = geostatsp::RFsimulate(trueParamIso,mydat)$sim

mydat$Y = -3 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
		mydat$U + rnorm(length(mydat), 0, sd=sqrt(trueParamAniso["nugget"]))

mydat$Yiso = -3 + 0.5*mydat$cov1 + 0.2*mydat$cov2 + 
    mydat$Uiso + rnorm(length(mydat), 0, sd=sqrt(trueParamIso["nugget"]))

mydat$Ybc = (mydat$Y*0.5+1)^2

mydat$YbcIso = (mydat$Yiso*0.5+1)^2

print(range(mydat$Ybc, na.rm=TRUE))
print(range(mydat$YbcIso, na.rm=TRUE))

date()
myres = likfitLgm(
    formula=Ybc ~ cov1 + cov2, 
    data=mydat,
    coordinates=mydat,
    param=c(range=0.01,nugget=2,shape=2, 
        anisoAngleDegrees=20, anisoRatio=2,
        boxcox=1), 
    paramToEstimate = c("range","nugget",
        "anisoRatio","anisoAngleDegrees",
        "boxcox","shape"),
    reml=TRUE
)
date()
myres$opt$logL
myres$par

date()
myresIso = likfitLgm(
    formula=YbcIso ~ cov1 + cov2, 
    data=mydat,
    param=c(shape=2,boxcox=0.4),
    paramToEstimate = c(
        "range","nugget",
        "boxcox","shape"),
    reml=TRUE
)
date()
myresIso$opt$logL
myresIso$par

myres$summary[,grep("^ci", colnames(myres$summary),invert=TRUE)]

loglikLgm(formula=Ybc ~ cov1 + cov2, 
    data=mydat, 
    param=myres$param
)

# only estimate variance
myres = likfitLgm(
    formula=Ybc ~ cov1 + cov2, 
    data=mydat,
    coordinates=mydat,
    param=c(range=0.01,nugget=2,shape=2, 
        anisoAngleDegrees=20, anisoRatio=2,
        boxcox=1), 
    paramToEstimate = c("variance"),
    reml=TRUE
)


pdf("ligfitLgm.pdf")
par(mfrow=c(1,2))

myraster = rast(nrows=30,ncols=30,xmin=0,xmax=1,ymin=0,ymax=1)
covEst = matern(myraster, y=c(0.5, 0.5), par=myres$param)
covTrue = matern(myraster, y=c(0.5, 0.5), par=trueParamAniso)

plot(covEst, main="estimate")
plot(covTrue, main="true")

dev.off()

if(interactive()  | Sys.info()['user'] =='patrick') {
  
  
  
library("geostatsp")
data("swissRain")
swissAltitude = unwrap(swissAltitude)

sr2 = unwrap(swissRain)
sr2$elev = extract(swissAltitude, sr2)
swissFitAgain = likfitLgm(data=sr2, 
		formula=rain~ elev,
		param=c(range=1000,shape=1,nugget=0.1,boxcox=0.5),
		paramToEstimate = c("range","nugget")
)
swissFitAgain$par		
}

