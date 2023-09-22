library("geostatsp")

param = c(range=1, shape=1.5,	anisoRatio=2, anisoAngleDegrees=-25)

matern(c(0, 0.001, 100000), param=param)

#x=c(0, 0.001, 100000);param=c(param, variance=1)

#resultFull = .C("matern", as.double(x), as.integer(length(x)),
#		as.double(param["range"]), as.double(param["shape"]),
#		as.double(param["variance"]))


# example with raster
myraster = rast(nrows=40,ncols=60,xmin=-3,xmax=3,ymin=-2,ymax=2)

# plot correlation of each cell with the origin
myMatern = matern(myraster, y=c(0,0), param=param)
myMatern[1:3,1:3]
 


bob = function(x) {
thepar = attributes(x)$param
pdf(tempfile("matern", tmpdir=".", fileext=".pdf"))
plot(x, main=
				paste(
				paste(names(thepar), thepar, sep="="),
				collapse=", "),cex.main=0.5
)
dev.off()
}

bob(myMatern)


bob(matern(myraster, y=c(1,-0.5), 
				param =	c(range=1, shape=1.5,	anisoRatio=2, anisoAngleDegrees=25)			
						)
)

bob(matern(myraster,y= c(0,0), 
				param =	c(range=1, shape=25.1,	anisoRatio=2, anisoAngleDegrees=-25)			
		)
)

bob(matern(myraster,y= c(0,0), 
				param =	c(range=0, shape=1.5,	anisoRatio=2, anisoAngleDegrees=-25)			
		)
)

bob(matern(myraster, y=c(0,0), 
				param =	c(range=100000, shape=1.5,	anisoRatio=2, anisoAngleDegrees=-25)			
		)
)


				
				

# correlation matrix for all cells with each other
myraster = rast(nrows=4,ncols=6,xmin=-3,xmax=3,ymin=-2,ymax=2)
myMatern = matern(myraster, param=c(range=0, shape=2))

dim(myMatern)
myMatern[1:3,1:3]

param = c(range=0.2, shape=1.5)
set.seed(0)

mypoints = vect(cbind(runif(10), runif(10)), 
		crs = "epsg:2000",
		atts=data.frame(id=1:10))

myDist = forceSymmetric(as.matrix(distance(mypoints)))

myDist3 = myDist[lower.tri(myDist)]

myMatern1 = matern(mypoints, param)
myMatern2 = matern(myDist, param)
myMatern3 = matrix(NA, length(mypoints),length(mypoints))
myMatern3[lower.tri(myMatern3)] = matern(myDist3, param)



class(myMatern1)
class(myMatern2)
class(myMatern3)


(myMatern1 - myMatern2)[1:8, 1:4]
sum((myMatern1 - myMatern2)^2)
(myMatern1 - myMatern3)[1:8, 1:4]
sum((myMatern1 - myMatern3)^2, na.rm=TRUE)


